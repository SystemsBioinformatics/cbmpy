'''
Created on Nov 11, 2014

@author: arne
'''

from __future__ import division, print_function
from __future__ import absolute_import

from math import isnan

from pyscescbm.solver import CBSolver
from .decomposition import Leaf
from .sparserationals import Matrix
from . import matroid

class EFMEnumerator():
    '''
    classdocs
    '''


    def __init__(self, mnet, node, excluded):
        """Creates a new EFMEnumerator for the given network
        
        We assume that all reactions in mnet are irreversible, i.e.
        that all lower flux bounds are greater or equal zero.
        """
        self.mnet = mnet
        self.node = node
        self.excluded = excluded
        self.module = node.getLeaves(excluded)
        # check if all reactions are in the module irreversible
        for r in self.module:
            assert mnet.getReactionLowerBound(r) >= 0
        
        # build stoichiometric matrix
        self.matrix = Matrix()
        labels, _ = self.matrix.addMetabolicNetwork(mnet)
        labidx = [i for i in range(len(labels)) if labels[i] in self.module]
        self.matcols = [labels[i] for i in labidx]
        self.matrix = self.matrix[:,labidx] # we only need the submatrix
        self.matroid = matroid.fromMatrix(self.matrix, self.matcols)
        
        self.tol = 1e-5
        self.minflux = 1e-4
        self.__buildLP()
        
    def __buildLP(self):
        """ build the LP for the computations """
        self.lp = CBSolver.createSolver(self.mnet)
        self.fixed = set()
        for r in self.mnet.reactions:
            # store which reactions have a fixed flux rate, because
            # those get a special treatment for minimality
            span = r.getFVAdata()[3]
            if (span < self.tol):
                self.fixed.add(r.getId())
#            if r.getId() not in self.module:
#                # remove all bounds from the variable
#                # unless it has a fixed flux rate
#                span = r.getFVAdata()[3]
#                if (span > self.tol or isnan(span)):
#                    # remove flux bounds
#                    self.lp.setLowerBounds({r.getId():None})
#                    self.lp.setUpperBounds({r.getId():None})
        # clear objective function
        self.lp.setObjective(reset=True)
        
#     def setDecomposition(self, decomposition):
#         """ Sets the decomposition to use for enumeration"""
#         self.decomp = decomposition
    
    def isFeasible(self, face):
        """ checks if face is a feasible face of the module 
        
        A face is feasible if there exists a point in the polyhedron 
        (maybe, we generalize this a bit?) where the 
        reactions of the face are the only reactions used from the reactions in
        the module.
        """
        bounds = {}
        for r in self.module:
            if r in face:
                # set lower bound to self.minflux
                lb = max(self.minflux, self.mnet.getReactionLowerBound(r))
                ub = self.mnet.getReactionUpperBound(r)
                if lb > ub:
                    return False
                bounds[r] = (lb, ub)
            else:
                # fix flux to zero
                lb = self.mnet.getReactionLowerBound(r)
                if lb > 0:
                    return False
                bounds[r] = (lb, 0)
 
        self.lp.setBounds(bounds)
        self.lp.solve()
        if self.lp.getSolutionStatus() != 'LPS_UNDEF':
            return self.lp.isFeasible()
        else:
            # try to resolve
            self.__buildLP()
            self.lp.setBounds(bounds)
            if self.lp.getSolutionStatus() != 'LPS_UNDEF':
                return self.lp.isFeasible()
            else:
                print('warning: unable to solve LP')
                return True # if we keep it, we don't miss solutions
            
        
    def isMinimal(self, face):
        """ checks if face is a minimal face of the module """
        #var = self.getVariable(face)
        var = face.difference(self.fixed) # this assumes that the face is feasible
        return self.matroid.isIndependent(var)
#         faceidx = [i for i in range(len(self.matcols)) if self.matcols[i] in var]
#         submatrix = self.matrix[:, faceidx]
#         m = matroid.fromMatrix(submatrix, faceidx)
#         if len(m.nbasis) == 0:
#             return True
#         else:
#             return False
    
    def getVariable(self, face):
        """run some kind of FVA on the reactions in the face
        
        we remark that reactions not in the face are fixed to 0 and hence 
        have no variability.
        """
        sol = self.lp.getSolution() # fetch old solution
        
        bounds = {}
        solfeasible = True
        for r in self.module:
            if r in face:
                # set lower bound to self.minflux
                lb = max(0, self.mnet.getReactionLowerBound(r))
                ub = self.mnet.getReactionUpperBound(r)
                bounds[r] = (lb, ub)
                # sol should never violate this
            else:
                # fix flux to zero
                lb = self.mnet.getReactionLowerBound(r)
                bounds[r] = (lb, 0)
                if sol[r] > self.tol:
                    solfeasible = False
        
        self.lp.setBounds(bounds)
        
        # if sol is still feasible, use it for init
        if solfeasible:
            maxFlux = sol
            minFlux = sol.copy()
        else:
            self.lp.solve()
            sol = self.lp.getSolution()
            maxFlux = sol
            minFlux = sol.copy()

        for r in face:
            if maxFlux[r] - minFlux[r] < self.tol:
                # solve min and max problem
                # do max first, since most lb are 0
                self.lp.setObjective(coef={r:1}, sense='max', reset=False)
                self.lp.solve()
                if not self.lp.isFeasible():
                    # this should not happen, but numerical instability etc.
                    # can cause this.
                    # Try to solve again
                    self.__buildLP()
                    self.lp.setBounds(bounds)
                    self.lp.setObjective(coef={r:1}, sense='max', reset=False)
                    self.lp.solve()
                    if not self.lp.isFeasible():
                        print('warning: unable to solve FVA lp!')
                sol = self.lp.getSolution()
                for s, v in sol.items():
                    if maxFlux[s] < v:
                        maxFlux[s] = v
                    if minFlux[s] > v:
                        minFlux[s] = v
                        
                # still not variable?
                if maxFlux[r] - minFlux[r] < self.tol:
                    # solve minimization
                    self.lp.setObjective(coef={r:1}, sense='min', reset=False)
                    self.lp.solve()
                    if not self.lp.isFeasible():
                        # this should not happen, but numerical instability etc.
                        # can cause this.
                        # Try to solve again
                        self.__buildLP()
                        self.lp.setBounds(bounds)
                        self.lp.setObjective(coef={r:1}, sense='min', reset=False)
                        self.lp.solve()
                        if not self.lp.isFeasible():
                            print('warning: unable to solve FVA lp!')
                    sol = self.lp.getSolution()
                    for s, v in sol.items():
                        if maxFlux[s] < v:
                            maxFlux[s] = v
                        if minFlux[s] > v:
                            minFlux[s] = v
                # fast reset objective
                self.lp.setObjective(coef={r:0}, sense='min', reset=False)
        
        variable = set()
        for r in face:
            if maxFlux[r] - minFlux[r] > self.tol:
                variable.add(r)
        return variable
    
    def enumerateVertices(self):
        """enumerate the vertices of the subnetwork reachable from node
        
        The subnetwork reachable from node consists of the leaves that can be
        reached in the decomposition from node without passing through excluded.
        """
       
    def test(self, face):
        face = face.intersection(self.module)
        assert self.isFeasible(face)
        assert self.isMinimal(face), '%s' % face.__str__() 
        
        adj = set(self.node.getAdjacent())
        for v in adj:
            if v != self.excluded:
                enum = EFMEnumerator(self.mnet, v, self.node)
                enum.test(face)
                
        return True
     
    def enumerateMinimal(self):
        """ enumerates minimal A-faces 
        
        otherwise the same as enumerateVertices.
        
        Note: In contrast to the paper, we represent faces by the reactions
        that are used (not blocked)
        """
        #module = node.getLeaves(excluded)
        if isinstance(self.node, Leaf):
            assert len(self.module) == 1
            F = set()
            if self.isFeasible(set()):
                F.add(frozenset())
            face = frozenset(self.module)
            if self.isFeasible(face):
                if self.isMinimal(face):
                    F.add(face)
            return F
        else:
            adj = set(self.node.getAdjacent())
            ###
            # for the moment lets not use any runtime optimizations
            ###
            v1 = adj.pop()
            while v1 == self.excluded:
                v1 = adj.pop()
            enum1 = EFMEnumerator(self.mnet, v1, self.node)
            F1 = enum1.enumerateMinimal()
            print('%s: %d' % (v1.__str__(), len(F1)) )
            
            # if we have multiple elements, just add them in arbitrary order
            # this corresponds to a caterpillar decomposition
            # when we merge the results we have to be careful what the base set
            # of elements of the module is, since we create virtual modules
            # in between.
            fullmodule = set(self.module)
            self.module = set(enum1.module)
            while len(adj) > 0:
                v2 = adj.pop()
                if v2 != self.excluded:
                    enum2 = EFMEnumerator(self.mnet, v2, self.node)
                    F2 = enum2.enumerateMinimal()
                    print('%s: %d' % (v2.__str__(), len(F2)) )
                    self.module.update(enum2.module)
                    F = set() # result set
                    print('combine %d * %d' % (len(F1), len(F2)) )
                    i = 0
                    for f1 in F1:
                        for f2 in F2:
                            f = f1.union(f2)
                            if self.isFeasible(f):
                                if self.isMinimal(f):
                                    F.add(f)
                        i = i+1
                        print('%d of %d (%d per step)' % (i, len(F1), len(F2)))
                    print('merge: %d' % len(F) )
                    F1 = F
            # at the end we should have reached the original module
            assert self.module == fullmodule
            return F1
            
            ###
            # potential, better code
            ###
#             assert self.excluded in adj
#             adj.discard(self.excluded)
#             if len(adj) > 2:
#                 # deal with parallel or coparallel elements
#                 if self.node.matroid.isParallel():
#                     # just enumerate for each and take the union
#                     F = set()
#                     for v in adj:
#                         F1 = self.enumerateMinimal(v, self.node)
#                         F.update(F1)
#                     return F
#                 else:
#                     assert self.node.matroid.isCoparallel()
#                     # either no reaction is used, or from every submodule
#                     # a reaction is used
#                     hasEmpty = False
#                     emptyset = frozenset()
#                     F = {emptyset}
#                     for v in adj:
#                         F1 = self.enumerateMinimal(v, self.node)
#                         if emptyset in F1:
#                             hasEmpty = True
#                             F1.discard(emptyset)
#                         Fnew = set()
#                         for f in F:
#                             for f1 in F1:
#                                 Fnew.add(f + f1)
#                         F = Fnew
#                         
#                     if hasEmpty:
#                         F.add(emptyset)
#                     return F
#             else:
#                 v1 = adj.pop()
#                 v2 = adj.pop()
#                 enum1 = EFMEnumerator(self.mnet, v1, self.node)
#                 enum2 = EFMEnumerator(self.mnet, v2, self.node)
#                 F1 = enum1.enumerateMinimal()
#                 F2 = enum2.enumerateMinimal()
#                 F = set() # result set
#                 for f1 in F1:
#                     for f2 in F2:
#                         f = f1 + f2
#                         if self.isFeasible(f):
#                             if self.isMinimal(f):
#                                 F.add(f)
#                 return F
                    