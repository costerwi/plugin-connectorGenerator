"""Abaqus CAE plugin to find and create couplings and wire features connecting
circular edges matching the size and Parts of the two edges selected by the user.

This is intended to make it much more convenient to fasten layers of midsurface
together using connectors when several circular fastener holes are available but
there are no fasteners.

Carl Osterwisch, October 2024
"""

from __future__ import print_function
from abaqus import *
from abaqusConstants import *
from bisect import bisect_left
import numpy as np
import os

DEBUG = os.environ.get('DEBUG')

def uniqueKey(repository, baseName='Item'):
    "Return a new unique key within the repository"
    n = 0
    while 0 == n or name in repository:
        n -= 1
        name = baseName + str(n)
    return name


def edgeId(edge):
    "Return a hashable identity for this edge"
    return edge.index, edge.instanceName


def referencePair(rp1, rp2):
    "Return unique identifier for these two refrence points"
    return tuple(sorted([rp1, rp2]))


def tangentEdges(edge, radius):
    "Return edgeArray of adjacent circular tangent edges of the specified radius"
    assert radius > 0
    tangentEdges = edge.getEdgesByEdgeAngle(1.0) # tangent edges within 1 degree
    adjacentTangents = tangentEdges[0:0]
    def addAdjacentTangents(newEdge):
        nonlocal adjacentTangents
        i = tangentEdges.index(newEdge)
        adjacentTangents += tangentEdges[i:i+1]
        for adjacentEdge in newEdge.getAdjacentEdges():
            if adjacentEdge in adjacentTangents:
                continue # already added
            if not adjacentEdge in tangentEdges:
                continue # not tangent
            try:
                if abs(radius - adjacentEdge.getRadius())/radius > 0.01:
                    continue # wrong radius
            except:
                continue # non-circular edge
            # TODO check for common center
            addAdjacentTangents(adjacentEdge)
    addAdjacentTangents(edge) # start from original edge
    return adjacentTangents


def getSimlarEdges(rootAssembly, edge0, radius):
    "Return list of edge arrays with same radius as edge0 in all instances of this Part"
    instance0 = rootAssembly.instances[edge0.instanceName]
    instance0Edges = [tangentEdges(edge0, radius)] # edge0 forms the first item in the list
    for edge in instance0.edges:  # search edges within original instance
        if edge.isReferenceRep:
            continue # reference geom
        if any(edge in edgeArray for edgeArray in instance0Edges):
            continue # already in the list
        try:
            if abs((radius - edge.getRadius())/radius) > 0.01:
                continue # wrong radius
            instance0Edges.append(tangentEdges(edge, radius))
        except:
            continue # ignore all errors

    # Add the same edges arrays from all instance of this part
    allSimilarEdges = instance0Edges.copy()
    for otherInstance in rootAssembly.instances.values():
        if not hasattr(otherInstance, 'partName'):
            continue # assembly instance
        if otherInstance.partName != instance0.partName:
            continue # different part
        if otherInstance.name == instance0.name:
            continue # same instance
        if rootAssembly.features[otherInstance.name].isSuppressed():
            continue # suppressed instance
        allSimilarEdges.extend(otherInstance.edges.getSequenceFromMask(
                edgeArray.getMask()) for edgeArray in instance0Edges)

    return allSimilarEdges


couplingPoints = {} # dict edgeId -> coupled rp
def reloadCouplings(model):
    "Rebuild couplingPoints based on couplings between reference points and edges"
    couplingPoints.clear()
    rootAssembly = model.rootAssembly
    rps = set(repr(rp) for rp in rootAssembly.referencePoints.values())
    remove = []
    for constraint in model.constraints.values():
        if not hasattr(constraint, 'couplingType'):
            continue # not a coupling
        if constraint.surface[1] != 'Assembly':
            continue
        surfaceName = constraint.surface[0]
        try:
            surface = rootAssembly.allInternalSurfaces[surfaceName]
        except KeyError:
            surface = rootAssembly.surfaces[surfaceName]
        if not surface.edges:
            continue # must have edges
        edge0 = min(surface.edges)
        controlSetName = constraint.controlPoint[0]
        try:
            controlSet = rootAssembly.allInternalSets[controlSetName]
        except KeyError:
            controlSet = rootAssembly.sets[controlSetName]
        if len(controlSet.referencePoints) != 1:
            continue # must have one reference point
        rp = controlSet.referencePoints[0]
        if repr(rp) in rps:
            couplingPoints[edgeId(edge0)] = rp
        else: # rp no longer exists
            remove.append(constraint.name)
    for constraintName in remove:
        print('Deleting constraint', constraintName, 'due to missing RP')
        del model.constraints[constraintName]


barePoints = {} # dict edgeId -> uncoupled rpFeature
def deleteUnusedCenters(rootAssembly):
    "Remove any new but uncoupled reference points"
    rootAssembly.deleteFeatures([feature.name for feature in barePoints.values()])
    barePoints.clear()


def centerPoint(model, edgeArray):
    "Return reference point at center of edgeArray"
    rootAssembly = model.rootAssembly
    # Search all existing couplings for this edge, return its referencePoint if found
    edge0 = min(edgeArray)
    edgeId0 = edgeId(edge0)
    if edgeId0 in couplingPoints:
        return couplingPoints[edgeId0]
    if edgeId0 in barePoints:
        rpFeature = barePoints[edgeId0]
        return rootAssembly.referencePoints[rpFeature.id]
    # Create a new reference point
    instance = rootAssembly.instances[edge0.instanceName]
    rpFeature = rootAssembly.ReferencePoint(point=instance.InterestingPoint(edge0, CENTER))
    barePoints[edgeId0] = rpFeature
    return rootAssembly.referencePoints[rpFeature.id]


def makeSpider(model, edgeArray):
    "Create coupling between edgeArray and its center rp"
    import regionToolset
    edge0 = min(edgeArray)
    edgeId0 = edgeId(edge0)
    if edgeId0 in couplingPoints:
        return # already connected by coupling
    rp = centerPoint(model, edgeArray)
    controlRegion = regionToolset.Region( referencePoints=[rp] )
    surfaceRegion = regionToolset.Region(side1Edges=edgeArray)
    coupling = model.Coupling(name=uniqueKey(model.constraints, 'CenterCoupling'),
        controlPoint=controlRegion,
        surface=surfaceRegion,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        rotationalCouplingType=ROTATIONAL_STRUCTURAL)
    couplingPoints[edgeId0] = rp
    del barePoints[edgeId0]
    return coupling


connectedPoints = []  # sorted list of connected point pairs

def reloadConnectedPoints(rootAssembly):
    "Refresh list connectedPoints based on existing model assembly features"
    connectedPoints.clear()
    groupByChild = {}
    for rpid, rp in rootAssembly.referencePoints.items():
        rpfeat = rootAssembly.featuresById[rpid]
        groupByChild.setdefault(rpfeat.children, []).append(rp)
    for rpList in groupByChild.values():
        if 2 != len(rpList):
            continue # not a 2 RP feature
        pair = referencePair(*rpList)
        connectedPoints.insert(bisect_left(connectedPoints, pair), pair) # insert sorted


def wireBetweenCenters(model, rpA, rpB):
    "Create a wire feature between center of edgeA and edgeB"

    rootAssembly = model.rootAssembly
    pair = referencePair(rpA, rpB)
    insertion = bisect_left(connectedPoints, pair)
    if insertion < len(connectedPoints) and pair == connectedPoints[insertion]:
        return None # wire already exists

    wire = rootAssembly.WirePolyLine(points=((rpA, rpB), ), meshable=False)
    connectedPoints.insert(insertion, pair) # insert sorted
    newName = uniqueKey(rootAssembly.features, 'WireConnector')
    rootAssembly.features.changeKey(fromName=wire.name, toName=newName)
    return rootAssembly.features[newName]


def addConnectors(edge1, edge2):
    "Main method called by CAE"

    from scipy.spatial import KDTree
    viewport = session.viewports[session.currentViewportName]
    rootAssembly = viewport.displayedObject
    model = mdb.models[rootAssembly.modelName]

    # Check for bad input from user
    edges = edge1, edge2
    radii = [edge.getRadius() for edge in edges] # will raise exception if not a radius
    if edge2 in tangentEdges(edge1, radii[0]):
        raise ValueError('The same edge was selected twice')
    partNames = [rootAssembly.instances[edge.instanceName].partName for edge in edges]

    reloadCouplings(model)
    deleteUnusedCenters(rootAssembly)
    reloadConnectedPoints(rootAssembly)

    try:
        viewport.disableColorCodeUpdates() # suspend updates for better performance

        similarEdges1 = getSimlarEdges(rootAssembly, edge1, radii[0])
        rp1 = [centerPoint(model, edgeArray) for edgeArray in similarEdges1]
        coords1 = [rootAssembly.getCoordinates(rp) for rp in rp1]

        # TODO handle special case of edge1.instanceName == edge2.instanceName => instances must always match

        if DEBUG:
            print('Parts', *partNames)
            print('Radii', *radii)
        if partNames[0] != partNames[1] or abs(radii[1] - radii[0])/radii[1] > 0.01:
            # edge1 and edge2 have different radii, different similarEdges
            if DEBUG:
                print('Non-matching parts or radii')
            similarEdges2 = getSimlarEdges(rootAssembly, edge2, radii[1])
            rp2 = [centerPoint(model, edgeArray) for edgeArray in similarEdges2]
            coords2 = [rootAssembly.getCoordinates(rp) for rp in rp2]
            boundDistance = np.linalg.norm(np.asarray(coords1[0]) - coords2[0]) + 0.05*sum(radii)
            # Find edge centers in similarEdges2 closest to centers of similarEdges1
            pointTree = KDTree(coords2)
            distances, index2 = pointTree.query(coords1, distance_upper_bound=boundDistance)
        else:
            # same part and radius for both edges; similarEdges2 will be same as similarEdges1
            if DEBUG:
                print('Matching parts and radii')
            similarEdges2 = similarEdges1
            rp2 = rp1
            coords2 = coords1
            for i, edgeList in enumerate(similarEdges2):
                if edge2 in edgeList:
                    boundDistance = np.linalg.norm(np.asarray(coords1[0]) - coords2[i]) + 0.1*radii[0]
                    break
            else:
                raise RuntimeError('Matching edge not found')
            # Find edge centers closest to each other
            pointTree = KDTree(coords2)
            distances, index2 = pointTree.query(
                    coords1,
                    k=[2],  # use second closest point to skip point matching with itself
                    distance_upper_bound=boundDistance)
            distances = distances.flatten()
            index2 = index2.flatten()

        row2Points = set() # keep track to prevent multiple edge1 edges connecting to the same edge2 edge
        wires = []
        for row in np.argsort(distances):
            if distances[row] > boundDistance:
                break # edge arrays in this row and remaining rows are missing matches
            row2 = index2[row]
            if row2 in row2Points:
                continue # another wire already connected to this edgeArray; shortest distance wins
            if rp1[row] == rp2[row2]:
                if DEBUG:
                    print('Skipping', row, row2)
                continue # don't connect a point to itself
            row2Points.add(row2)
            makeSpider(model, similarEdges1[row1])
            makeSpider(model, similarEdges2[row2])
            try:
                wires.append(wireBetweenCenters(model, rp1[row], rp2[row2]))
            except Exception as e:
                print(repr(e))
                print(edgeId(similarEdges1[row]), edgeId(similarEdges2[row2]))
                print(repr(rp1[row]), repr(rp2[row2]))
        deleteUnusedCenters(rootAssembly)
        wireNames = set(w.name for w in wires if w is not None)
        newEdges = rootAssembly.edges[0:0]  # empty edgeArray
        for edge in rootAssembly.edges:
            if not edge.featureName in wireNames:
                continue
            newEdges += rootAssembly.edges[edge.index:edge.index + 1]

    finally:
        viewport.enableColorCodeUpdates() # enable viewport updates even if exception

    if not newEdges:
        print('No new wires added from {0[0]} edge diameter {1[0]:.3g} to {0[1]} edge diameter {1[1]:.3g}'.format(
            partNames, 2*np.asarray(radii)))
    else:
        name = uniqueKey(rootAssembly.sets, 'WireConnectors')
        rootAssembly.Set(name=name, edges=newEdges)
        print(len(newEdges), name,
            'added from {0[0]} edge diameter {1[0]:.3g} to {0[1]} edge diameter {1[1]:.3g}'.format(
            partNames, 2*np.asarray(radii)))
