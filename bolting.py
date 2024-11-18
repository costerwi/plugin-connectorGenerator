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
import regionToolset
import numpy as np
from scipy.spatial import KDTree

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
    "Return edgeArray of tangent edges with the same radius as edge"
    edgeArray = edge.getEdgesByEdgeAngle(1.0) # tangent edges within 1 degree
    for edge in edgeArray:
        if abs(radius - edge.getRadius())/radius > 0.01:
            raise ValueError('Tangent edges have inconsistent radii')
    return edgeArray


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


newPoints = []
edgeCenters = {} # dict edge.index -> rp
def resetCenters(model):
    "Rebuild edgeCenters based on couplings between reference points and edges"
    newPoints.clear()
    edgeCenters.clear()
    rootAssembly = model.rootAssembly
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
        edgeCenters[edgeId(edge0)] = controlSet.referencePoints[0]


def deleteUnusedCenters(rootAssembly):
    "Remove the new but unusued reference points"
    unusedCenterPoints = []
    for featureId in newPoints:
        feature = rootAssembly.featuresById[featureId]
        if not feature.children: # not used by a coupling
            unusedCenterPoints.append(feature.name)
    rootAssembly.deleteFeatures(unusedCenterPoints)


def centerPoint(model, edgeArray):
    "Create reference point at center of circle and connect with coupling"
    rootAssembly = model.rootAssembly
    # Search all existing couplings for this edge, return its referencePoint if found
    edge0 = min(edgeArray)
    existing = edgeCenters.get(edgeId(edge0))
    if existing:
        return existing
    # Create a new reference point
    instance = rootAssembly.instances[edge0.instanceName]
    rpFeature = rootAssembly.ReferencePoint(point=instance.InterestingPoint(edge0, CENTER))
    newPoints.append(rpFeature.id)
    return rootAssembly.referencePoints[rpFeature.id]


def makeSpider(model, edgeArray, rp):
    "Create coupling between edgeArray and rp"
    edge0 = min(edgeArray)
    if edgeId(edge0) in edgeCenters:
        return # already connected by coupling
    controlRegion = regionToolset.Region( referencePoints=[rp] )
    surfaceRegion = regionToolset.Region(side1Edges=edgeArray)
    coupling = model.Coupling(name=uniqueKey(model.constraints, 'CenterCoupling'),
        controlPoint=controlRegion,
        surface=surfaceRegion,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        rotationalCouplingType=ROTATIONAL_STRUCTURAL)
    edgeCenters[edgeId(edge0)] = rp # remember
    return coupling


connectedPoints = []  # sorted list of existing connector point pairs

def resetConnectedPoints(rootAssembly):
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
    newName = uniqueKey(rootAssembly.features, 'WireBolt')
    rootAssembly.features.changeKey(fromName=wire.name, toName=newName)
    return rootAssembly.features[newName]


def addConnectors(edge1, edge2):
    "Main method called by CAE"

    viewport = session.viewports[session.currentViewportName]
    rootAssembly = viewport.displayedObject
    model = mdb.models[rootAssembly.modelName]

    # Check for bad input from user
    if edge1.instanceName == edge2.instanceName:
        raise ValueError('Edges are from the same instance')
    edges = edge1, edge2
    radii = [edge.getRadius() for edge in edges] # will raise exception if not a radius
    for edge in edges:
        instance = rootAssembly.instances[edge.instanceName]
        assert hasattr(instance, 'partName')

    resetCenters(model)
    resetConnectedPoints(rootAssembly)

    similarEdges1 = getSimlarEdges(rootAssembly, edge1, radii[0])
    similarEdges2 = getSimlarEdges(rootAssembly, edge2, radii[1])

    try:
        viewport.disableColorCodeUpdates() # suspend updates for better performance

        rp1 = [centerPoint(model, edgeArray) for edgeArray in similarEdges1]
        rp2 = [centerPoint(model, edgeArray) for edgeArray in similarEdges2]

        coords1 = [rootAssembly.getCoordinates(rp) for rp in rp1]
        coords2 = [rootAssembly.getCoordinates(rp) for rp in rp2]

        boundDistance = np.linalg.norm(np.asarray(coords1[0]) - coords2[0]) + 0.5*sum(radii)

        # Find edge centers in similarEdges2 closest to centers of similarEdges1
        pointTree = KDTree(coords2)

        distances, index2 = pointTree.query(coords1, distance_upper_bound=boundDistance)

        row2Points = set() # keep track to prevent multiple edge1 edges connecting to the same edge2 edge
        wires = []
        for row in np.argsort(distances):
            if distances[row] > boundDistance:
                break # edge arrays in this row and remaining rows are missing matches
            row2 = index2[row]
            if row2 in row2Points:
                continue # another wire already connected to this edgeArray; shortest distance wins
            row2Points.add(row2)
            makeSpider(model, similarEdges1[row], rp1[row])
            makeSpider(model, similarEdges2[row2], rp2[row2])
            wires.append(wireBetweenCenters(model, rp1[row], rp2[row2]))
        deleteUnusedCenters(rootAssembly)
        wireNames = set(w.name for w in wires if w is not None)
        newEdges = rootAssembly.edges[0:0]  # empty edgeArray
        for edge in rootAssembly.edges:
            if not edge.featureName in wireNames:
                continue
            newEdges += rootAssembly.edges[edge.index:edge.index + 1]

    finally:
        viewport.enableColorCodeUpdates() # enable viewport updates even if exception

    partNames = [rootAssembly.instances[edge.instanceName].partName for edge in edges]
    if not newEdges:
        print('No new wires added from {0[0]} hole diameter {1[0]:.3g} to {0[1]} hole diameter {1[1]:.3g}'.format(
            partNames, 2*np.asarray(radii)))
    else:
        name = uniqueKey(rootAssembly.sets, 'BoltWires')
        rootAssembly.Set(name=name, edges=newEdges)

        print(len(newEdges), name,
            'added from {0[0]} hole diameter {1[0]:.3g} to {0[1]} hole diameter {1[1]:.3g}'.format(
            partNames, 2*np.asarray(radii)))
