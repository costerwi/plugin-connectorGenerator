from abaqus import *
from abaqusConstants import *
import regionToolset
import numpy as np
from scipy.spatial import KDTree


def uniqueKey(base, repository):
    n = 0
    while 0 == n or name in repository:
        n -= 1
        name = base + str(n)
    return name


def getSimlarEdges(rootAssembly, edge0):
    "Find similar size circular edges in all instances of this Part"
    radius = edge0.getRadius()
    similarEdges = [edge0] # Similar size edges within the original instance
    instance = rootAssembly.instances[edge0.instanceName]
    for edge in instance.edges:
        if edge.isReferenceRep:
            continue # reference geom
        if edge.index == edge0.index:
            continue # original edge
        if len(edge.getVertices()) > 1:
            continue # incomplete circle TODO consider circles made up of multiple edges
        #f = edge.getFaces()
        #if len(f) > 1 or f[0] != faces[0]:
        #    continue # in case of sketch editing
        try:
            r2 = edge.getRadius()
        except Exception as E:
            continue # not recognized as an arc
        if abs((radius - r2)/radius) > 0.01:
            continue # wrong radius
        similarEdges.append(edge)

    # Add the same edges from all instance of this part
    copiedEdges = similarEdges.copy()
    for otherInstance in rootAssembly.instances.values():
        if not hasattr(otherInstance, 'partName'):
            continue # assembly instance
        if rootAssembly.features[otherInstance.name].isSuppressed():
            continue # suppressed instance
        if otherInstance.name == instance.name:
            continue # same instance
        if otherInstance.partName != instance.partName:
            continue # different part
        similarEdges.extend(otherInstance.edges[e.index] for e in edges)

    return copiedEdges


def hashEdge(edge):
    "Return unique tuple for this edge"
    return edge.instanceName, edge.index


def centerPoint(model, edge):
    "Create reference point at center of circle and connect with coupling"
    rootAssembly = model.rootAssembly
    # Search all existing couplings for this edge, return its referencePoint
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
        if len(surface.edges) != 1:
            continue
        if hashEdge(edge) != hashEdge(surface.edges[0]):
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
    instance = rootAssembly.instances[edge.instanceName]
    rpFeature = rootAssembly.ReferencePoint(point=instance.InterestingPoint(edge, CENTER))
    rp = rootAssembly.referencePoints[rpFeature.id]
    rootAssembly.features.changeKey(fromName=rpFeature.name,
                                    toName=uniqueKey('BoltRP', rootAssembly.features))
    controlRegion = regionToolset.Region( referencePoints=[rp] )
    surfaceRegion = regionToolset.Region(side1Edges=instance.edges[edge.index:edge.index+1])
    #surfaceRegion = regionToolset.Region(side1Edges=([edge],) )
    #surfaceRegion = regionToolset.Region(side1Edges=[instance.edges[edge.index]] )
    model.Coupling(name=uniqueKey('BoltCoupling', model.constraints)),
        controlPoint=controlRegion,
        surface=surfaceRegion,
        influenceRadius=WHOLE_SURFACE,
        couplingType=KINEMATIC,
        rotationalCouplingType=ROTATIONAL_STRUCTURAL)
    return rp


newConnections = set()  # make sure an edge is only added once
def bond(model, edgeA, edgeB):
    "Create a wire feature between center of edgeA and edgeB"
    for edge in edgeA, edgeB:
        if hashEdge(edge) in newConnections:
            return None # Avoid duplicate connections
    rootAssembly = model.rootAssembly
    rpA = centerPoint(model, edgeA)
    rpB = centerPoint(model, edgeB)
    for edge in edgeA, edgeB:
        newConnections.add( hashEdge(edge) )
    wire = rootAssembly.WirePolyLine(points=((rpA, rpB), ), mergeType=IMPRINT, meshable=False)
    newName = uniqueKey('Bolt-{}-{}'.format(edgeA.instanceName, edgeB.instanceName), rootAssembly.features)
    rootAssembly.features.changeKey(fromName=wire.name, toName=newName)
    return rootAssembly.features[newName]


def myMethod(edge1, edge2):
    # Check for bad input from user
    viewport = session.viewports[session.currentViewportName]
    rootAssembly = viewport.displayedObject
    model = mdb.models[rootAssembly.modelName]

    if edge1.instanceName == edge2.instanceName:
        raise ValueError('Edges are from the same instance')
    edges = edge1, edge2
    radius = [edge.getRadius() for edge in (edge1, edge2)] # will raise exception if not a radius
    if any(r < 1e-3 for r in radius):
        raise ValueError('Radius is very small')
    for edge in edges:
        if len(edge.getVertices()) > 1:
            raise ValueError('Edge {} forms incomplete circle'.format(edge.instanceName))
        instance = rootAssembly.instances[edge.instanceName]
        assert hasattr(instance, 'partName')

    newConnections.clear() # reset

    wires = [bond(model, edge1, edge2)]

    similarEdges1 = getSimlarEdges(rootAssembly, edge1)
    similarEdges2 = getSimlarEdges(rootAssembly, edge2)

    # Find and connect similar sized circular edges on all instances of the two parts
    pointTree = KDTree([e.pointOn[0] for e in similarEdges2])
    distances, index2 = pointTree.query([e.pointOn[0] for e in similarEdges1])
    for edge, distance, j in zip(similarEdges1, distances, index2):
        if distance > 20: # TODO base this limit on the distance between points in the initial bond
            continue
        wires.append(bond(model, edge, similarEdges2[j]))

    wireNames = set(w.name for w in wires if w is not None)
    geomsequence = None
    for index, edge in enumerate(rootAssembly.edges):
        if not edge.featureName in wireNames:
            continue
        if geomsequence is None:
            geomsequence = rootAssembly.edges[index:index + 1]
        else:
            geomsequence += rootAssembly.edges[index:index + 1]
    boltSet = rootAssembly.Set(name=uniqueKey('BoltWires', rootAssembly.sets),
                     edges=geomsequence)

    print(len(geomsequence), boltSet.name, 'wires D{:.3g} to D{:.3g}'.format(2*radius[0], 2*radius[1]))
