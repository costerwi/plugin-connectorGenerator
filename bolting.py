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


def tangentEdges(edge):
    "Return edgeArray of tangent edges with the same radius as edge"
    radius = edge.getRadius()
    edgeArray = edge.getEdgesByEdgeAngle(1.0) # tangent edges within 1 degree
    for edge in edgeArray:
        r = edge.getRadius()
        if abs(r - radius)/radius > 0.01:
            raise ValueError('Tangent edges have inconsistent radii')
    return edgeArray


def getSimlarEdges(rootAssembly, edge0):
    "Return list of edge arrays with same radius as edge0 in all instances of this Part"
    radius = edge0.getRadius()
    instance0 = rootAssembly.instances[edge0.instanceName]
    instance0Edges = [tangentEdges(edge0)] # edge0 forms the first item in the list
    for edge in instance0.edges:  # search edges within original instance
        if edge.isReferenceRep:
            continue # reference geom
        if any(edge in edgeArray for edgeArray in instance0Edges):
            continue # already in the list
        try:
            r2 = edge.getRadius()
        except:
            continue # not recognized as an arc
        if abs((radius - r2)/radius) > 0.01:
            continue # wrong radius
        try:
            instance0Edges.append(tangentEdges(edge))
        except:
            continue # failed to construct tangent edges

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


def hashEdge(edge):
    "Return unique tuple for this edge"
    return edge.instanceName, edge.index


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


connectedPoints = []  # make sure an edge is only added once
def wireBetweenEdgeCenters(model, edgeArrayA, edgeArrayB):
    "Create a wire feature between center of edgeA and edgeB"

    rootAssembly = model.rootAssembly
    rpA = centerPoint(model, edgeArrayA)
    rpB = centerPoint(model, edgeArrayB)
    pairId = referencePair(rpA, rpB)
    if pairId in connectedPoints:
        return None # wire already exists

    wire = rootAssembly.WirePolyLine(points=((rpA, rpB), ), meshable=False)
    connectedPoints.append(pairId)
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
        connectedPoints.append(referencePair(*rp))

    similarEdges1 = getSimlarEdges(rootAssembly, edge1)
    similarEdges2 = getSimlarEdges(rootAssembly, edge2)

    distance0 = np.linalg.norm(np.asarray(edge1.pointOn[0]) - edge2.pointOn[0])
    maxDistance = 1.5*(distance0 + sum(radii))

    # Find edges in similarEdges2 closest to simpiarEdges1
    pointTree = KDTree([min(edgeArray).pointOn[0] for edgeArray in similarEdges2])
    distances, index2 = pointTree.query([min(edgeArray).pointOn[0] for edgeArray in similarEdges1],
                                        distance_upper_bound=maxDistance)

    row2Points = set() # keep track to prevent multiple wires connecting to the same edge2 similar edge
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
    geomsequence = None
    for index, edge in enumerate(rootAssembly.edges):
        if not edge.featureName in wireNames:
            continue
        if geomsequence is None:
            geomsequence = rootAssembly.edges[index:index + 1]
        else:
            geomsequence += rootAssembly.edges[index:index + 1]

    if not geomsequence:
        print('No new wires added from hole diameter {:.3g} to hole diameter {:.3g}'.format(
            2*radii[0], 2*radii[1]))
    else:
        name = uniqueKey(rootAssembly.sets, 'BoltWires')
        rootAssembly.Set(name=name, edges=geomsequence)

        print(len(geomsequence), name,
            'added from hole diameter {:.3g} to hole diameter {:.3g}'.format(
            2*radii[0], 2*radii[1]))
