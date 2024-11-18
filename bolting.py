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


def centerPoint(model, edgeArray):
    "Create reference point at center of circle and connect with coupling"
    rootAssembly = model.rootAssembly
    # Search all existing couplings for this edge, return its referencePoint if found
    edge0 = min(edgeArray)
    for constraint in model.constraints.values():
        if not hasattr(constraint, 'couplingType'):
            continue # not a coupling
        if constraint.surface[1] != 'Assembly':
            continue
        surfaceName = constraint.surface[0]
        if surfaceName in rootAssembly.allInternalSurfaces:
            surface = rootAssembly.allInternalSurfaces[surfaceName]
        else:
            surface = rootAssembly.surfaces[surfaceName]
        if not edge0 in surface.edges:
            continue # not the same edge
        controlSetName = constraint.controlPoint[0]
        if controlSetName in rootAssembly.allInternalSets:
            controlSet = rootAssembly.allInternalSets[controlSetName]
        else:
            controlSet = rootAssembly.sets[controlSetName]
        if len(controlSet.referencePoints) < 1:
            continue # no reference points
        # Matching! Return the existing controlPoint
        return controlSet.referencePoints[0]

    # Create a new reference point and coupling constraint
    instance = rootAssembly.instances[edge0.instanceName]
    rpFeature = rootAssembly.ReferencePoint(point=instance.InterestingPoint(edge0, CENTER))
    rp = rootAssembly.referencePoints[rpFeature.id]
    rootAssembly.features.changeKey(fromName=rpFeature.name,
                                    toName=uniqueKey(rootAssembly.features, 'BoltRP'))
    controlRegion = regionToolset.Region( referencePoints=[rp] )
    surfaceRegion = regionToolset.Region(side1Edges=edgeArray)
    model.Coupling(name=uniqueKey(model.constraints, 'BoltCoupling'),
        controlPoint=controlRegion,
        surface=surfaceRegion,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        rotationalCouplingType=ROTATIONAL_STRUCTURAL)
    return rp


connectedPoints = []  # sorted list of existing connector point pairs
def wireBetweenEdgeCenters(model, edgeArrayA, edgeArrayB):
    "Create a wire feature between center of edgeA and edgeB"

    rootAssembly = model.rootAssembly
    rpA = centerPoint(model, edgeArrayA)
    rpB = centerPoint(model, edgeArrayB)
    pair = referencePair(rpA, rpB)
    insertion = bisect_left(connectedPoints, pair)
    if insertion < len(connectedPoints) and pair == connectedPoints[insertion]:
        return None # wire already exists

    wire = rootAssembly.WirePolyLine(points=((rpA, rpB), ), meshable=False)
    connectedPoints.insert(insertion, pair) # insert sorted
    newName = uniqueKey(rootAssembly.features,
            'Bolt-{}-{}'.format(min(edgeArrayA).instanceName, min(edgeArrayB).instanceName))
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

    # Collect connectedPoints of existing wires to avoid duplication
    connectedPoints.clear()
    for edge in rootAssembly.edges:
        vertexList = edge.getVertices()
        if 2 != len(vertexList):
            continue # not a 2-point connector
        vertices = (rootAssembly.vertices[i] for i in vertexList)
        rp = (rootAssembly.referencePoints.findAt(*v.pointOn) for v in vertices) # TODO is findAt most efficient?
        pair = referencePair(*rp)
        connectedPoints.insert(bisect_left(connectedPoints, pair), pair) # insert sorted

    similarEdges1 = getSimlarEdges(rootAssembly, edge1, radii[0])
    similarEdges2 = getSimlarEdges(rootAssembly, edge2, radii[1])

    distance0 = np.linalg.norm(np.asarray(edge1.pointOn[0]) - edge2.pointOn[0])
    maxDistance = 1.5*(distance0 + sum(radii))

    # Find edges in similarEdges2 closest to simpiarEdges1
    pointTree = KDTree([min(edgeArray).pointOn[0] for edgeArray in similarEdges2])
    distances, index2 = pointTree.query([min(edgeArray).pointOn[0] for edgeArray in similarEdges1],
                                        distance_upper_bound=maxDistance)

    row2Points = set() # keep track to prevent multiple edge1 edges connecting to the same edge2 edge
    wires = []
    for row in np.argsort(distances):
        if distances[row] > maxDistance:
            break # edge arrays in this row and remaining rows are too far away from each other
        row2 = index2[row]
        if row2 in row2Points:
            continue # already connected to this edgeArray; shortest distance wins
        row2Points.add(row2)
        wires.append(wireBetweenEdgeCenters(model, similarEdges1[row], similarEdges2[row2]))

    wireNames = set(w.name for w in wires if w is not None)
    newEdges = rootAssembly.edges[0:0]  # empty edgeArray
    for edge in rootAssembly.edges:
        if not edge.featureName in wireNames:
            continue
        newEdges += rootAssembly.edges[edge.index:edge.index + 1]


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
