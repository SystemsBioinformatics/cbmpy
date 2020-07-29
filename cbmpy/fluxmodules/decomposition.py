"""
CBMPy: fluxmodules decomposiiton module
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Arne C. Reimers
Contact email: arne.c.reimers@gmail.com
"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import numpy

from .matroid import fromMatrix
from .random_color import generate_new_color

class Edge:
    """ Every edge connects two vertices """
    
    def __init__(self, a, b):
        """ initialize edge connecting vertices a and b """
        self.a = a
        self.b = b
        
    def getOther(self, vertex):
        """ give the other vertex then the specified one 
        
        The specified vertex must be one of this edges vertices
        """
        if self.a == vertex:
            return self.b
        else:
            assert self.b == vertex, 'vertex must be a vertex of this edge'
            return self.a

class AbstractVertex:
    """ A vertex in a graph. It can be a leaf or an inner node """

    maxId = 0
    """ for giving unique ids for hashing """

    def __init__(self):
        self.__id = AbstractVertex.maxId
        AbstractVertex.maxId += 1
        pass
        
    def id(self):
        return self.__id
    
    def __hash__(self):
        return self.__id
    
    def __eq__(self, other):
        """ two abstract vertices are only the same if they are identical """
        return self is other
    
    def label(self):
        print("subclasses must override this")

    def replaceAdjacent(self, old, new):
        """ replace the adjacent node old by new """
        print('subclasses must override this')
        assert False

    def getAdjacent(self):
        """ list adjacent nodes """
        print('subclasses must override this')
        assert False

    def getLeaves(self, fromVertex):
        """ enumerates all the leaves of this vertex that are reachable without
        using fromVertex
        
        fromVertex must be a vertex adjacent to this vertex
        """
        print('subclasses must override this')
        assert False

class Vertex(AbstractVertex):
    """ represents an inner vertex (i.e., not a leaf)
    
    Every Vertex contains a map (edges) of adjacent edges resp. vertices.
    Each Edge has a list of elements associated to it (the interface).
    Each Vertex stores a matroid, whoose elements are contributed by the edges.
    
    The map edges has the following form:
    The keys are the incident vertices.
    The values are the list of elements contributed by the edge
    """
    def __init__(self):
        AbstractVertex.__init__(self)
        self.edges={}

    def __repr__(self):
        return "N%s" % self.id()
    
    def __str__(self):
        return "N%s" % self.id()
    
    def label(self):
        return self.__str__()

    def setUpRoot(self, matroid):
        """ set ups this vertex as the root node for the given matroid.
        
        This means we add leaves for all the elements of the matroid to this
        vertex
        """
        for e in matroid.elems:
            l = Leaf(e, self)
            self.edges[l] = [e]
        self.matroid = matroid

    def splitOff(self, adjacent):
        """ splits of a list of adjacent vertices.
        
        adjacent must be a list of vertices that are adjacent to this node.
        The ist adjacent must contain more than one element 
        """
        assert len(adjacent) >= 2, 'splitting off would have no effect'
        
        # the next assert may be useful, but it can make things also much
        # more complicated and it does not really hurt the branch decompositon
        #assert len(adjacent) <= len(self.edges) <= 2
        
        # collect elements to split off
        toSplit = set([])
        for a in adjacent:
            toSplit.update(self.edges[a])
        # find a basis of self.matroid such that the only basis elements of
        # toSplit are used by non-basis elements of toSplit
        
        # in the following loop, we modify the set of non-basis elements,
        # so first copy it to avoid concurrent modification errors
        for e in list(self.matroid.nbasis):
            if e not in toSplit:
                c = self.matroid.fundamentalCircuit(e)
                intersection = c.intersection(toSplit)
                if len(intersection) >= 1:
                    # put e into the basis 
                    # and one element of the intersection into nbasis
                    # note that the new non-basis element is in toSplit, so we
                    # don't have to do anything with it
                    self.matroid.exchange(e,intersection.pop())
       
        # our matroid is now in the desired form:
        # no fundamental circuit of a nbasis element outside toSplit uses basis
        # elements of toSplit
        # and no cofundamental circuit  of a basis element in toSplit uses nbasis
        # elements outside of toSplit

        # create new vertex
        vsplit = Vertex()
        for a in adjacent:
            # move edge from self to vsplit
            vsplit.edges[a] = self.edges[a]
            # delete adjacency of a
            del self.edges[a]
            a.replaceAdjacent(self, vsplit)
            
        # the dual matroid is of the desired form for the complement of toSplit,
        # which we will use to simplify the interface of vsplit
        vsplit.matroid = self.matroid.dual() 
        
        # reduce representation of toSplit for self
        interface = self.__simplifyInterface(toSplit)
        self.edges[vsplit] = interface # add vsplit as neighbour
        
        if len(adjacent) >= 3 or True:
            # reduce representation of the complement of toSplit for vsplit
            co_toSplit = set(vsplit.matroid.elems) - toSplit
            cointerface = vsplit.__simplifyInterface(co_toSplit)
            vsplit.edges[self] = cointerface # add self as neighbour to vsplit
            # in the end we actually want to have the primal matroid stored
            vsplit.matroid = vsplit.matroid.dual()
        else:
            # vsplit is a node with 3 vertices, we will probably not touch it
            # again, so be lazy and do not calculate the cointerface
            vsplit.edges[self] = len(interface) * [None]
            vsplit.matroid = None
        
        # return the new vertex, maybe we want to do something with it
        return vsplit
        
    def replaceAdjacent(self, old, new):
        assert old in self.edges, 'old must be an adjacent vertex'
        # copy interface
        self.edges[new] = self.edges[old]
        # delete old edge
        del self.edges[old]
        
     
    def __simplifyInterface(self, toSplit):
        """ reduces the interface of the elements to split off
        
        For this method we assume:
        No fundamental circuit of a non-basis element outside of toSplit must
        contain a basis element inside toSplit
        """
        assert self.__checkSplit(toSplit), 'basis has wrong form'
        
        self.matroid = self.matroid.contract(
                        toSplit.intersection(self.matroid.basis)
                                            )
        # this gives us an interface of the split
        interface = []
        for i in self.matroid.nbasis:
            if i in toSplit:
                interface.append(i)
        # however, we want to simplify this interface
        # it should only have as many elements as the connectivity is
        
        # take the reduced representation matrix of the non-basis elements in 
        # toSplit in self.matroid and interpret it as a direct representation 
        # matrix and find a basis
        toDelete = [i for i in self.matroid.nbasis if i not in toSplit]
        # we only delete non-basis elements, so delete does not mess with basis
        submat = self.matroid.delete(toDelete) 
        mred = fromMatrix(submat.rep, interface)
        # the basis elements of mred represent the interface
        # so we delete all the other elements
        toDelete = [i for i in interface if i in mred.nbasis]
        self.matroid = self.matroid.delete(toDelete)
        
        return mred.basis
        
    def __checkSplit(self, toSplit):
        for e in list(self.matroid.nbasis):
            if e not in toSplit:
                c = self.matroid.fundamentalCircuit(e)
                if not c.isdisjoint(toSplit):
                    return False
        return True

    def getAdjacent(self):
        return frozenset(self.edges.keys())

    def getLeaves(self, fromVertex):
        """ enumerates all the leaves of this vertex that are reachable without
        using fromVertex
        
        If fromVertex is not a vertex incident to this vertex, then all leaves
        are enumerated.
        """
        leaves = set()
        for v in self.edges.keys():
            if v != fromVertex:
                leaves.update(v.getLeaves(self))
        return leaves
    
    def isSimple(self):
        if len(self.matroid.elems) == 0:
            return True
        # check if its parallel elements
        parallel = self.matroid.parallel()
        if len(parallel) == 1:
            par = parallel.pop()
            if len(par) == len(self.matroid.elems):
                assert par == frozenset(self.matroid.elems)
                return True
                
        # check if its cpparallel elements
        coparallel = self.matroid.coparallel()
        if len(coparallel) == 1:
            copar = coparallel.pop()
            if len(copar) == len(self.matroid.elems):
                assert copar == frozenset(self.matroid.elems)
                return True
            
        return False

class Leaf(AbstractVertex):

    def __init__(self, element, adjacent):
        """ creates leaf representing an matroid element
        
        element is the matroid element to represent
        adjacent is the vertex to which this leaf is added
        """
        AbstractVertex.__init__(self)
        self.elem = element
        self.adj = adjacent

    def __repr__(self):
        return "Leaf%s <%s>" % (self.id(), self.elem)

    def __str__(self):
        return "Leaf%s <%s>" % (self.id(), self.elem)
    
    def label(self):
        return self.elem.__str__()

    def replaceAdjacent(self, old, new):
        assert self.adj == old, 'old must be single adjacent node of this leaf'
        self.adj = new

    def getAdjacent(self):
        return frozenset([self.adj])

    def getLeaves(self, fromVertex):
        assert self.adj == fromVertex, 'fromEdge must be from the leaf'
        return {self.elem}

class Decomposition:
    """ Hierarchical branch decomposition 
    
    Actually can also represent a partial branch decomposition as intermediate
    results.
    
    It has the following fields:
    matroid - the matroid of this branch decomposition
    root - a vertex of the graph representing the branch decomposition
    
    TODO: more details
    """
    
    def __init__(self, matroid):
        """ initializes as star branch decomposition """
        self.matroid = matroid
        self.root = Vertex()
        self.root.setUpRoot(matroid)

    def listNonLeaves(self):
        """ list all non-leaf vertices """
        return self.listNonLeavesRec(self.root, None)

    def listNonLeavesRec(self, vertex, source):
        if isinstance(vertex, Leaf):
            return set()
        
        result = {vertex}
        for v in vertex.getAdjacent():
            if v != source:
                result.update(self.listNonLeavesRec(v, vertex))
        return result
        
    def getVertex(self, vertexid):
        """ returns the vertex with the given id
        
        If no vertex with the given id exists, None is returned.
        """
        for v in self.listNonLeaves():
            if v.id() == vertexid:
                return v
        return None

    def splitZeroModules(self, vertex):
        """ branch given vertex into 0 modules 
        
        Requires that all adjacent vertices have at most 1-dimensional interface
        """
        lookUpAdjacent = {}
        for adjacent, elems in vertex.edges.items():
            assert len(elems) <= 1
            for e in elems:
                lookUpAdjacent[e] = adjacent
                
        modules = vertex.matroid.findModules()
        if len(modules) >= 2: # otherwise, splitting helps not
            for m in modules:
                if len(m.elems) > 1: # don't do anything for loops or coloops
                    toSplit = set([])
                    for e in m.elems:
                        toSplit.add(lookUpAdjacent[e])
                        
                    vertex.splitOff(toSplit)
                else:
                    #for loops and coloops we have to fix the interface
                    e = m.elems[0]
                    v = lookUpAdjacent[e]
                    vertex.edges[v] = []
                    vertex.matroid = vertex.matroid.delete([e])
                    
    def splitSimple(self):
        """ split all vertices into 0-modules and 1-modules """
        for v in self.listNonLeaves():
            self.splitZeroModules(v)
        
        for v in self.listNonLeaves():
            self.splitOneModules(v)
                        
    def splitOneModules(self, vertex):
        """ branch given vertex into 1 modules 
        
        TODO:we have to deal with grouped elements
        """
        # first try to group parallel and coparallel elements
        didSomething = True
        while didSomething:
            didSomething = False
            # create lookup table for adjacent vertices
            # the matroid changes, so we update it for each iteration
            lookUpAdjacent = {}
            for adjacent, elems in vertex.edges.items():
                for e in elems:
                    lookUpAdjacent[e] = adjacent
                
            # find parallel elements
            parallel = vertex.matroid.parallel()
            print("%s parallel: %s" % (vertex, parallel))
            for group in parallel:
                assert len(group) >= 2, 'groups must be non-trivial'
                # we split vertices not elements, so get the vertices
                toSplit = set()
                for e in group:
                    v = lookUpAdjacent[e]
                    if len(vertex.edges[v]) <= 1:
                        toSplit.add(v)
                if (len(toSplit) >= 2) and (len(toSplit) <= len(vertex.edges)-2):
                    vertex.splitOff(toSplit)
                    didSomething = True
            
            # update lookup table if we changed something
            if didSomething:
                lookUpAdjacent = {}
                for adjacent, elems in vertex.edges.items():
                    for e in elems:
                        lookUpAdjacent[e] = adjacent

            # find coparallel elements
            coparallel = vertex.matroid.coparallel()
            print("%s coparallel: %s" % (vertex, coparallel))
            for group in coparallel:
                assert len(group) >= 2, 'groups must be non-trivial'
                # we split vertices not elements, so get the vertices
                toSplit = set()
                for e in group:
                    v = lookUpAdjacent[e]
                    if len(vertex.edges[v]) <= 1:
                        toSplit.add(v)
                if (len(toSplit) >= 2) and (len(toSplit) <= len(vertex.edges)-2):
                    vertex.splitOff(toSplit)
                    didSomething = True
        
        print('TODO: more complicated 1-modules')
    
    def makeSubcubic(self):
        """ split all nodes that have degree >= 4 in an arbitrary manner."""
        for v in self.listNonLeaves():
            assert isinstance(v, Vertex)
            while len(v.getAdjacent()) >= 4:
                # split off two arbitrary nodes
                adj = v.getAdjacent()
                aiter = adj.__iter__()
                v.splitOff([aiter.next(), aiter.next()])
                
    
    def isParallelSeries(self):
        """ checks if all nodes have only parallel or only coparallel elems.
        
        This does not check if the matroid is parallel-series. 
        """
        for v in self.listNonLeaves():
#             ok = False
#             # check if its parallel elements
#             parallel = v.matroid.parallel()
#             if len(parallel) == 1:
#                 par = parallel.pop()
#                 if len(par) == len(v.matroid.elems):
#                     assert par == frozenset(v.matroid.elems)
#                     ok = True
#                     
#             # check if its cpparallel elements
#             coparallel = v.matroid.coparallel()
#             if len(coparallel) == 1:
#                 copar = coparallel.pop()
#                 if len(copar) == len(v.matroid.elems):
#                     assert copar == frozenset(v.matroid.elems)
#                     ok = True
#             
#             if not ok:
            if not v.isSimple():
                return False
            
        # all nodes are either parallel or coparallel
        return True
    
    def verifyFullyBranched(self):
        """ verifies that this branch decomposition is fully branched
        
        The branch decomposition is considered fully branched if all vertices
        satisfy:
        - is degree 3 or
        - all elements are parallel
        - all elements are coparallel
        
        If successful, True is returned.
        Otherwise, a NotFullyBranchedException is raised
        """
        for v in self.listNonLeaves():
            ok = False
            if len(v.edges) <= 3:
                ok = True
            
            if not ok:
                if v.isSimple():
                    ok = True
#                 # check if its parallel elements
#                 parallel = v.matroid.parallel()
#                 if len(parallel) == 1:
#                     par = parallel.pop()
#                     if len(par) == len(v.matroid.elems):
#                         assert par == frozenset(v.matroid.elems)
#                         ok = True
#                     
#             if not ok:
#                 # check if its cpparallel elements
#                 coparallel = v.matroid.coparallel()
#                 if len(coparallel) == 1:
#                     copar = coparallel.pop()
#                     if len(copar) == len(v.matroid.elems):
#                         assert copar == frozenset(v.matroid.elems)
#                         ok = True
            
            if not ok:
                raise NotFullyBranchedException(v)
                return False
            
        # all nodes are either parallel or coparallel
        return True
        
    def verifyEdgeWidth(self):
        """ verifies if the edge with really coincides with the interfaces.
        
        In case of success, true is returned.
        Otherwise an EdgeWidthException is raised.
        """
        for v in self.listNonLeaves():
            for (w, interface) in v.edges.items():
                sep = w.getLeaves(v)
                print(sep)
                conn = self.matroid.connectivity(sep)
                if len(interface) != conn:
                    raise EdgeWidthException(v,w,conn)
                    return False
        return True
    
    def getWidth(self):
        """ computes the width of this branch decomposition.
        """
        width = 0
        for v in self.listNonLeaves():
            for (w, interface) in v.edges.items():
                width = max(width, len(interface))
        return width
    
    def poolmanMethod(self, vertex, similarity=None):
        """ decomposes the matroid using the method suggested by Poolman et al.
        
        This method iteratively groups these elements together that are most 
        similar. Similarity by sets of nodes is done via the mean.
        For more details see doi:10.1016/j.jtbi.2007.08.005
        
        This method takes a similarity matrix as input, which allows it to 
        also apply the method to other similarity matrices then those suggested
        in the original work. 
        If no similarity matrix is given, then a row similarity matrix is used,
        which corresponds to the matrix suggested by Poolman et al.
        """

        if similarity == None:
            similarity = self.matroid.similarityMatrix("rows")
        
                
        # make a copy so that we don't make a mess with the input
        sim = similarity.matrix.copy()
        
        if set(similarity.elems) == set(vertex.matroid.elems):
            # similarity matrix is directly for the vertex matroid
            
            # merge columns representing vertices
            # preallocate, so we don't get indexing problems
            vertexIndices = len(similarity.elems) * [None]
            for w,elems in vertex.edges.items():
                # typically len(elems) == 1
                sumrow = numpy.zeros(len(similarity.elems))
                i = -1
                for e in elems:
                    i = similarity.elems.index(e)
                    sumrow += sim[i,:]
#                    sumrow = numpy.maximum(sumrow, sim[i,:])
                    sim[i,:] = -1 # does not overwrite important entries
                    sim[:,i] = -1 # does not overwrite important entries
                    sumrow[i] = -1 # make sure the value is negative
                
                if len(elems) >= 1:
                    assert i >= 0
                    sumrow /= len(elems)
                    sim[i,:] = sumrow
                    sim[:,i] = sumrow
                    vertexIndices[i] = (w, len(elems))
        else:
            assert set(similarity.elems) == set(self.matroid.elems)
            # similarity is of the whole-network matroid
            
            # merge elements in adjacent nodes
            # preallocate, so we don't get indexing problems
            vertexIndices = len(similarity.elems) * [None]
            adj = vertex.getAdjacent()
            for a in adj:
                tomerge = a.getLeaves(vertex)
                sumrow = numpy.zeros(len(similarity.elems))
                i = -1
                for e in tomerge:
                    i = similarity.elems.index(e)
                    sumrow += sim[i,:]
#                    sumrow = numpy.maximum(sumrow, sim[i,:])
                    sim[i,:] = -1 # does not overwrite important entries
                    sim[:,i] = -1 # does not overwrite important entries
                    sumrow[i] = -1 # make sure the value is negative
                
                if len(tomerge) >= 1:
                    assert i >= 0
                    sumrow /= len(tomerge)
                    sim[i,:] = sumrow
                    sim[:,i] = sumrow
                    vertexIndices[i] = (a, len(tomerge))
            
        
        # now we can start with the actual algorithm
        while len(vertex.edges) >= 4:
            (i,j) = numpy.unravel_index(sim.argmax(),sim.shape)
            assert i != j
            assert sim[i,j] >= 0
            assert vertexIndices[i] != None
            assert vertexIndices[j] != None
            (vi, ni) = vertexIndices[i]
            (vj, nj) = vertexIndices[j]
            # merge i and j, because they are most similar
            vsplit = vertex.splitOff([vi, vj])
            # create entry for vertex in similarity matrix
            newrow = (ni*sim[i,:] + nj*sim[j,:])/(ni+nj)
#            newrow = numpy.maximum(sim[j,:], sim[i,:])
            sim[i,:] = newrow
            sim[:,i] = newrow
            sim[j,:] = -1
            sim[:,j] = -1
            sim[i,i] = -1
            vertexIndices[i] = (vsplit, ni + nj)
            vertexIndices[j] = None
            
    def makeGraphViz(self, cmod=None):
        """ create graphviz string to draw this decomposition """
        out = "graph G {\n"
        found = set()
        for v in self.listNonLeaves():
            for (w,e) in v.edges.items():
                if w not in found:
                    out += v.label() + " -- " + w.label()
                    out += " [label = "+len(e).__str__()+"]\n"
            found.add(v)
        
        # color leaves according to subsystems
        if cmod != None:
            subsystemColors = {} # to store which subsystem gets which color
            for w in self.matroid.elems:
                r = cmod.getReaction(w)
                if (r != None) and (r.hasAnnotation("SUBSYSTEM")):
                    subsystem = r.getAnnotation("SUBSYSTEM")
                    if not subsystemColors.has_key(subsystem):
                        # create new color
                        c = generate_new_color(subsystemColors.values())
                        subsystemColors[subsystem] = c
                    c = subsystemColors[subsystem]
                    out += w.__str__() + " [color=\"%f %f %f\"]\n" % c
            # make legend
            out += "subgraph cluster_legend {\n"
            out += "label=\"Legend\";\n"
            colorcount = 0
            for (subsystem, color) in subsystemColors.items():
                out += "subsystemcolor%d" % colorcount
                out += " [label = \"%s\" " % subsystem
                out += "color = \"%f %f %f\"]\n" % color
                colorcount += 1
            out += "}\n"
        out += "}\n"
        return out
    
    def write(self, filename):
        """ writes this decomposition to file """
        fo = open(filename, mode='w')
        width = self.getWidth()
        variables = len(self.matroid.elems)
        rank = len(self.matroid.basis)
        fo.write("vars %d\n" % variables)
        fo.write("rank %d\n" % rank)
        fo.write("width %d\n" % width)
        found = set()
        for v in self.listNonLeaves():
            for (w,e) in v.edges.items():
                if w not in found:
                    fo.write("%s %s %d\n" % (v.label(), w.label(), len(e)))
            found.add(v)
        fo.close()

class NotFullyBranchedException(Exception):
    """ Thrown if the branch decomposition is not fully branched."""
    def __init__(self, vertex):
        self.vertex = vertex
        
    def __str__(self):
        return "%s is not fully branched." % self.vertex
    
class EdgeWidthException(Exception):
    
    def __init__(self, vertex1, vertex2, conn):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.conn = conn
        
    def __str__(self):
        out = "connectivity %d\n" % self.conn
        if isinstance(self.vertex1, Vertex):
            e1 = self.vertex1.edges[self.vertex2]
            out += "%s to %s has: %s\n" % (self.vertex1, self.vertex2, e1)
        if isinstance(self.vertex2, Vertex):
            e1 = self.vertex1.edges[self.vertex2]
            out += "%s to %s has: %s\n" % (self.vertex1, self.vertex2, e1)
        return out