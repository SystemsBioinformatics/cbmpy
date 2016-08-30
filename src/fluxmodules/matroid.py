"""
CBMPy: fluxmodules matroid module
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

from math import ceil

import time
#import sympy
import numpy
#from sympy.core.cache import clear_cache

from .sparserationals import Matrix

class Matroid:
    """A linear matroid
    
    It has the fields:
    elems - lists the names of the elements
    basis - lists the names of the basis elements
    nbasis - lists the names of the non-basis elements (dual basis)
    rep - reduced representation matrix, 
            rows correspond to basis-elements and
            columns correspond to non-basis elements
    """
    #def __init__(self, matrix, elems):
    #    """create matroid from matrix storing elements as columns"""
    #    self.elems = __cols
    #    self.matrix = matrix
    #    self.buildRepresentation()
        
    def __init__(self, rep, basis, nbasis):
        """create matroid from given reduced representation matrix.
        
        The names of the rows are given by basis.
        The names of the columns are given by nbasis. 
        """
        self.elems = basis + nbasis
        self.basis = basis
        self.nbasis = nbasis
        self.rep = rep
    
            
    def findModules(self):
        """computes the separators (modules) of this matroid
        
        Each module is returned as a matroid.
        """
        usedBasis = set()
        usedNBasis = set()
        modules = []
        for r in range(len(self.basis)):
            if r not in usedBasis:
                c = self.__dfs(r)
                usedBasis.update(c['basis'])
                usedNBasis.update(c['nbasis'])
#                 m = Matroid(self.rep.extract(c['basis'], c['nbasis']),
#                             [self.basis[i] for i in c['basis']],
#                             [self.nbasis[i] for i in c['nbasis']]
#                            )
                m = Matroid(self.rep[c['basis'], c['nbasis']],
                            [self.basis[i] for i in c['basis']],
                            [self.nbasis[i] for i in c['nbasis']]
                           )
                modules.append(m)
        # by construction each basis element is contained in one submatroid
        # however, this need not be the case for the non-basis elements
        
        # each non-basis element that we did not find is a loop
        for r in range(len(self.nbasis)):
            if r not in usedNBasis:
                #m = Matroid(sympy.zeros(0,1), [], [r])
                m = Matroid(Matrix(), [], [r])
                modules.append(m)
        
        #clear_cache()
        
        return modules
        
    def __dfs(self, start):
        """performs a dfs to find the connected component of basis elem start"""
        basisFound = set([start])
        nbasisFound = set()
        todo = [start]
        while len(todo) > 0:
            r = todo.pop()
            for i in range(len(self.nbasis)):
                # if r is in the cycle and the cycle has not been added, 
                # then add all elements of the cycle to todo
                if (self.rep[r,i] != 0) and (i not in nbasisFound):
                    nbasisFound.add(i)
                    for s in range(len(self.basis)):
                        if (self.rep[s,i] != 0) and (s not in basisFound):
                            basisFound.add(s)
                            todo.append(s)
                            
        return {'basis':basisFound, 'nbasis':nbasisFound}
         
    def dual(self):
        """ returns the dual matroid 
        
        By this operation the original coefficient matrix is lost for the dual.
        """
        # we also care about the signs, so multiply by -1
        m = Matroid(self.rep.transpose(-1), list(self.nbasis), list(self.basis))
        return m
        
    def copy(self):
        """ copies this matroid """
        m = Matroid(self.rep.copy(), list(self.basis), list(self.nbasis))
        if hasattr(self, 'matrix'):
            m.matrix = self.matrix.copy()
        return m
       
    def contract(self, toContract):
        """ computes self / toContract
        
        The matroid with toContract contracted is returned. 
        This instance is not modified.
        """
        # put as many elements from toDelete into basis
        self.rank(toContract)
        retain_nbasis = [i for i in range(len(self.nbasis)) 
                            if self.nbasis[i] not in toContract]
        retain_basis = [i for i in range(len(self.basis)) 
                            if self.basis[i] not in toContract]
        basis = [self.basis[i] for i in retain_basis]
        nbasis = [self.nbasis[i] for i in retain_nbasis]
        
        submat = self.rep[retain_basis, retain_nbasis]
        m = Matroid(submat, basis, nbasis)
        #m = self.dual().delete(toContract).dual()
        if hasattr(self,'matrix'):
            print('TODO: do also the contraction on the matrix')
            
        #clear_cache()
        return m
        
    def exchange(self, a, b):
        """ exchanges a basis element with a non-basis element.
        
        requirement is that the fundamental circuit of the non-basis element
        uses the basis element
        """
        if a in self.basis:
            be = a
            assert b in self.nbasis, 'one element must be a non-basis element'
            nbe = b
        else:
            assert b in self.basis, 'one element must be a basis element'
            be = b
            assert a in self.nbasis, 'a, b must be elements of the matroid'
            nbe = a
        
        bidx = self.basis.index(be)
        nbidx = self.nbasis.index(nbe)
        nbidxcol = self.rep[:,nbidx]
        
        # We can simply switch be with nbe and only have to normalize the 
        # column. This is easy, so we do it last
        div = self.rep[bidx, nbidx]
        assert div != 0, 'basis elem and non-basis elem do not share circuit'
        
        # There may be more columns that use the basis element
        # since the basis element changed, we have to update
        for j in range(len(self.nbasis)):
            if (j != nbidx) and (self.rep[bidx,j] != 0):
                # add self.rep[:,nbidx] to cancel 
                # self.rep[bidx,j]
                coef = -self.rep[bidx,j] / div
                self.rep.coladd(j, nbidxcol, coef)
                #self.rep[:,j] += self.rep[:,nbidx]*coef
                assert self.rep[bidx,j] == 0
                # remember that we changed the meaning of the row,
                # so update the entry also
                self.rep[bidx,j] = coef
        
        # No do the switch with be and nbe
        self.rep[:,nbidx] /= div
        self.rep[bidx, nbidx] = 1/div
        
        # and also update the column and row names
        self.basis[bidx] = nbe
        self.nbasis[nbidx] = be
    
    def delete(self, toDelete):
        """ computes self \ toDelete
        
        The matroid with toDelete deleted is returned. 
        No elements of this instance are deleted. However, the current basis 
        can change.
        """
        
        # put as many elements from toDelete into non-basis
        self.corank(toDelete)
        retain_nbasis = [i for i in range(len(self.nbasis)) 
                            if self.nbasis[i] not in toDelete]
        retain_basis = [i for i in range(len(self.basis)) 
                            if self.basis[i] not in toDelete]
        basis = [self.basis[i] for i in retain_basis]
        nbasis = [self.nbasis[i] for i in retain_nbasis]
        
        submat = self.rep[retain_basis, retain_nbasis]
        m = Matroid(submat, basis, nbasis)
        
#         #######################################################################
#         ## old, borken code
#         ######################################################################
#         
#         # deletion of non-basis elements is easy
#         delete_nbasis = set([]) # columns to delete
#         for i in range(len(self.nbasis)):
#             if self.nbasis[i] in toDelete:
#                 delete_nbasis.add(i)
#                         
#         # for deletion of basis elements, we first have to make them non-basis 
#         # elements unless they are coloops
#         delete_basis = set([]) # for coloops
#         basis_new = list(self.basis) # new basis
#         for i in range(len(self.basis)):
#             if self.basis[i] in toDelete:
#                 switch_elem = -1
#                 for j in range(len(self.nbasis)):
#                     if (self.rep[i,j] != 0) and (j not in delete_nbasis):
#                         switch_elem = j
#                 if switch_elem == -1:
#                     # i is a coloop and we can just delete it
#                     delete_basis.add(i)
#                 else:
#                     basis_new[i] = self.nbasis[switch_elem]
#                     delete_nbasis.add(switch_elem)
#                     switch_col = self.rep[:,switch_elem]
#                     # now there may be more columns that use the basis element
#                     # since the basis element changed, we have to update
#                     # We do not use self.exchange(..) because we can save some
#                     # work.
#                     for j in range(len(self.nbasis)):
#                         if (self.rep[i,j] != 0) and (j not in delete_nbasis):
#                             # add self.rep[:,switch_elem] to cancel 
#                             # self.rep[i,j]
#                             coef = -self.rep[i,j] / switch_col[i,0]
#                             #coef = -self.rep[i,j] / self.rep[i,switch_elem]
#                             self.rep.coladd(j, switch_col, coef)
#                             #self.rep[:,j] += coef*self.rep[:,switch_elem]
#                             assert self.rep[i,j] == 0
#                             # remember that we changed the meaning of the row,
#                             # so update the entry also
#                             self.rep[i,j] = coef
#         
#         retain_nbasis = []
#         nbasis = []
#         for i in range(len(self.nbasis)):
#             if i not in delete_nbasis:
#                 retain_nbasis.append(i)
#                 nbasis.append(self.nbasis[i])
#                
#         retain_basis = []
#         basis = []
#         for i in range(len(self.basis)):
#             if i not in delete_basis:
#                 retain_basis.append(i)
#                 basis.append(basis_new[i])
#         
#         
#         #submat = self.rep.extract(retain_basis, retain_nbasis)
#         submat = self.rep[retain_basis, retain_nbasis]
#         m = Matroid(submat, basis, nbasis)
        
        # if we have the original matrix representation, also delete from that
        if hasattr(self, 'matrix'):
#             m.matrix = self.matrix.extract(
#                 range(self.matrix.rows),
#                 [i for i in range(len(self.elems)) 
#                     if self.elems[i] not in toDelete]
#                                           )
            m.matrix = self.matrix[
                range(self.matrix.rows()),
                [i for i in range(len(self.elems)) 
                    if self.elems[i] not in toDelete]
                                  ]
            # here the order of elements for sure gets mixed up, so set anew
            m.elems = [i for i in self.elems if i not in toDelete]
        
        #clear_cache()
        return m
    
    def parallel(self):
        """ list all groups of parallel elements 
        
        The result is a set of frozen sets, each set listing the ids of the 
        parallel elements
        """
        # basis elements are never parallel, unless there exists a
        # non-basis element with a single non-zero entry for the basis element
        # therefore, we just iterate through the non-basis elements
        par_groups = set([])
        
        found = set([])
        
        # to speed up the calculation, compute hash-values of the columns
        hashes = []
        for i in range(len(self.nbasis)):
            norm_factor = 0
            mul_factor = 73
            hashv = 1
            for j in range(len(self.basis)):
                v = self.rep[j,i] 
                if v != 0:
                    if norm_factor == 0:
                        norm_factor = 1/v
                    v = v*norm_factor
                    hashv += v.__hash__()
                hashv *= mul_factor
            hashes.append(hashv)
        
        for i in range(len(self.nbasis)):
            # only check if we did not already found it to be parallel to some
            # previous reaction
            if i not in found:
                par_group = [self.nbasis[i]]
                # check if parallel to some basis element
                basis_idx = -1
                basis_val = 0
                num_basis = 0
                for j in range(len(self.basis)):
                    if (num_basis < 2) and (self.rep[j,i] != 0):
                        num_basis += 1
                        if basis_idx == -1: # first non-zero entry
                            basis_idx = j 
                            basis_val = self.rep[j,i]
                        
                if num_basis == 0:
                    # its a loop and continue
                    continue
                if num_basis == 1:
                    par_group.append(self.basis[basis_idx])
                    
                assert basis_val != 0
                # check if parallel to other non-basis elements
                for j in range(len(self.nbasis))[i+1:]:
                    # if two cols are the same, their hashes must be also
                    if hashes[i] == hashes[j]:
                        coef = self.rep[basis_idx, j] / basis_val
                        # check if coef != 0 not only for performance reasons
                        # this also filters out loops
                        if coef != 0:
                            diff = self.rep[:,j]
                            diff.coladd(0, self.rep[:,i], -coef)
                            if diff.isZero():
                                par_group.append(self.nbasis[j])
                                found.add(j)
            
                # add par_group if non-trivial
                if len(par_group) >= 2:
                    par_groups.add(frozenset(par_group))
        
        return par_groups
        
    def coparallel(self):
        """ list all groups of coparallel elements 
        
        In metabolic network terminology the coparallel elements are the
        fully coupled reactions.
        
        The result is a set of frozensets, each set listing the ids of the 
        coparallel elements
        """
        return self.dual().parallel()
    
    def isParallel(self):
        parallel = self.parallel()
        if len(parallel) != 1:
            return False
        
        par = parallel.pop()
        if len(par) == len(self.matroid.elems):
            return True
        else:
            return False
    
    def isCoparallel(self):
        coparallel = self.coparallel()
        if len(coparallel) != 1:
            return False
        
        par = coparallel.pop()
        if len(par) == len(self.matroid.elems):
            return True
        else:
            return False
    
    def fundamentalCircuit(self, e):
        """ returns the fundamental circuit of an element
        
        If e is a non-basis element, the fundamental circuit, i.e. the
        support of the column of e, is returned.
        
        If e is a basis element, the cofundamental circuit, i.e. the support of
        the row of e, is returned.
        
        The returned fundamental circuit does not contain e itself
        """
        fcircuit = set([])
        
        if e in self.basis:
            # return cofundamental circuit (row)
            idx = self.basis.index(e)
            for j in range(len(self.nbasis)):
                if self.rep[idx,j] != 0:
                    fcircuit.add(self.nbasis[j])
        else:
            assert e in self.nbasis, 'this works only for elements of the matroid'
            # return fundamental circuit (column)
            idx = self.nbasis.index(e)
            for j in range(len(self.basis)):
                if self.rep[j,idx] != 0:
                    fcircuit.add(self.basis[j])
        return fcircuit
    
    def fundamentalCircuitValues(self, e):
        """ returns the fundamental circuit of an element with the coefficients
        
        If e is a non-basis element, the fundamental circuit, i.e. the
        support of the column of e, is returned.
        
        If e is a basis element, the cofundamental circuit, i.e. the support of
        the row of e, is returned.
        
        The returned fundamental circuit does not contain e itself. If needed,
        it should be added manually afterwards with a coefficient of 1.
        """
        fcircuit = {}
        
        if e in self.basis:
            # return cofundamental circuit (row)
            idx = self.basis.index(e)
            for j in range(len(self.nbasis)):
                if self.rep[idx,j] != 0:
                    fcircuit[self.nbasis[j]] = self.rep[idx,j]
        else:
            assert e in self.nbasis, 'this works only for elements of the matroid'
            # return fundamental circuit (column)
            idx = self.nbasis.index(e)
            for j in range(len(self.basis)):
                if self.rep[j,idx] != 0:
                    fcircuit[self.basis[j]] = self.rep[j,idx]
        return fcircuit
    
    def rank(self, elems):
        """ computes the rank of elems
        
        This method puts as many elements of elems into the basis. The number of
        elements from elems in the final basis is returned.
        """
        # try to make as many nbasis elements of elems basis-elements
        for i in range(len(self.nbasis)):
            if self.nbasis[i] in elems:
                # make it a basis element
                for j in range(len(self.basis)):
                    if (self.rep[j,i] != 0) and (self.basis[j] not in elems):
                        self.exchange(self.basis[j], self.nbasis[i])
                        break
        
        # count basis elements in elems
        numBasis = 0
        for i in range(len(self.basis)):
            if self.basis[i] in elems:
                numBasis += 1
        return numBasis

    def corank(self, elems):
        """ computes the corank of elems
        
        This method puts as many elements of elems into the non-basis. The 
        number of elements from elems in the final non-basis is returned.
        """
        # try to make as many nbasis elements of elems basis-elements
        for i in range(len(self.basis)):
            if self.basis[i] in elems:
                # make it a nbasis element
                for j in range(len(self.nbasis)):
                    if (self.rep[i,j] != 0) and (self.nbasis[j] not in elems):
                        self.exchange(self.nbasis[j], self.basis[i])
                        break
        
        # count non-basis elements in elems
        numNBasis = 0
        for i in range(len(self.nbasis)):
            if self.nbasis[i] in elems:
                numNBasis += 1
        return numNBasis

    
    def isIndependent(self, elems):
        """ checks if the given set of elements is independent.
        
        This does not modify this matroid, but may change its internal basis.
        If this method returns true, the current basis contains all elements of
        elems
        """
        # no independent set can have more elements than any basis
        if len(elems) > len(self.basis):
            return False 
        
        # we check by making all elements of elems basis elements
        for i in range(len(self.nbasis)):
            if self.nbasis[i] in elems:
                # make it a basis element
                found = False
                for j in range(len(self.basis)):
                    if (self.rep[j,i] != 0) and (self.basis[j] not in elems):
                        self.exchange(self.basis[j], self.nbasis[i])
                        found = True
                        break
                if not found:
                    return False
        # all elements of elems are now basis elements
        return True
    
    def isCoIndependent(self, elems):
        """ checks if the given set of elements is coindependent.
        
        This does not modify this matroid, but may change its internal basis.
        If this method returns true, the current basis contains all no elements 
        of elems
        
        This is the same as isIndependent, only on the dual
        """
        # no independent set can have more elements than any basis
        if len(elems) > len(self.nbasis):
            return False 
        
        # we check by making all elements of elems non-basis elements
        for j in range(len(self.basis)):
            if self.basis[j] in elems:
                # make it a basis element
                found = False
                for i in range(len(self.nbasis)):
                    if (self.rep[j,i] != 0) and (self.nbasis[i] not in elems):
                        self.exchange(self.basis[j], self.nbasis[i])
                        found = True
                        break
                if not found:
                    return False
        # all elements of elems are now non-basis elements
        return True
        
    def connectivity(self, separation):
        """ evaluate the connectivity function for the given separation.
        
            separation is a set of elements encoding one side of the separation.
            The other side is implied.
        """
        # for the computation we do the following:
        # put as many elements as possible from separation into the basis
        # this gives us the rank of separation.
        # We do the same for the complement.
        
        # make a copy, to reduce the number of pivots
        copy = self.copy()
        r1 = self.rank(separation)
        r2 = copy.rank(set(self.elems).difference(separation))
        r = len(self.basis)
        #clear_cache()
        return r1 + r2 - r
    
    def getInterface(self, separation):
        """ compute a list of circuits that define the interface of the
        given separation """
        
        # theoretically the following would be sufficient:
        #self.rank(separation)
        # but we do something more complicated...
        
        # Before we start the computation, we observe that eventually it is good 
        # to select circuits that are as simple as possible, so that the
        # selected interface vectors become as simple as possible.
        # typically the biomass reactions induce very_ugly circuits, hence we
        # do not want them to be in the basis, because then the fundamental
        # circuits get ugly.
        # HACK: We identify biomass reactions by their reaction NAME!!!
        # TODO: get some cleaner solution for the future!
        
        # find non-biomass reactions in the separation
        separation_non_bio = [r for r in separation if r.lower().find('bio') == -1]
        # find biomass reactions in the separation
        separation_bio = [r for r in separation if r.lower().find('bio') != -1]
        
        assert len(separation_non_bio) + len(separation_bio) == len(separation)
        
        self.rank(separation)
        
        # try to turn as many non-bio nbasis elements into basis elements by
        # exchanging with bio elements
        # use the property that every fundamental circuit of nbasis separation 
        # elements only contains separation elements
        for i in range(len(self.nbasis)):
            if self.nbasis[i] in separation and self.nbasis[i].find('bio') == -1:
                # make it a basis element
                for j in range(len(self.basis)):
                    if (self.rep[j,i] != 0) and (self.basis[j].lower().find('bio') != -1):
                        assert self.basis[j] in separation
                        self.exchange(self.basis[j], self.nbasis[i])
                        break
        
#        print(separation_bio)
#        self.rank(separation_bio)
#        print(self.basis)
#        print(separation_non_bio)
#        self.rank(separation_non_bio)
#        print(self.basis)
        # a maximal amount of separation elements are in the basis and of that
        # as few as possible are biomass reactions
        
        # Verify this!
        for r in self.nbasis:
            if r in separation:
                c = self.fundamentalCircuit(r)
                for s in c:
                    assert s in separation, (r,c,separation)
        
        # find biomass reactions not in the separation
        coseparation_bio = [r for r in self.basis if r not in separation and r.lower().find('bio') != -1]
        # ideally these reactions should not be in the basis, so try to
        # put a maximal amout of them into the nbasis
#        print(coseparation_bio)
        self.corank(coseparation_bio)
#        print(self.basis)
        # still the property is satisfied that a maximal amount of separation
        # elments are in the basis
        
        # we now have the following properties:
        # every fundamental circuit of an nbasis element in separation uses
        # only basis elements in separation
        # every cofundamental cocircuit of a basis element outside separation
        # uses only non-basis elements outside separation
        
        # Verify this!
        for r in self.nbasis:
            if r in separation:
                c = self.fundamentalCircuit(r)
                for s in c:
                    assert s in separation
        
        # we only have to deal with the fundamental circuits of nbasis elements
        # outside of separation.
        # They already define the interface, but the representation might be
        # redundant
        # Hence, we want to minimize it
                
        # Since the basis elements outside of the separation are of no interest,
        # we would simply contract them. However, this would make it more
        # difficult to reconstruct the original circuits.
        # There, we manually extract the reduced representation and take the
        # submatrix with rows of elements in separation and columns not in
        # separation
        rows = [i for i in range(len(self.basis)) if self.basis[i] in separation]
        cols = [i for i in range(len(self.nbasis)) if self.nbasis[i] not in separation]
        colNames = [self.nbasis[i] for i in cols]
        
        submat = self.rep[rows,cols]
        # redundant circuits are linearly dependent on the existing circuits
        # hence, we just need to find a basis of the reduced matrix
        # this we do by simply building the matroid
        mred = fromMatrix(submat, colNames)
        
        # here we continue our hack about the biomass reactions.
        # Theoretically, all bases of mred are equally good, but since we want
        # to have simple circuits, we prefer to have no biomass reactions in
        # the selected circuits and thus not in the basis of mred
#        print(mred.basis)
        col_bio = [r for r in colNames if r.lower().find('bio') != -1]
        mred.corank(col_bio)
#        print(mred.basis)
        
        interface = []
        for c in mred.basis:
            circuit = self.fundamentalCircuitValues(c)
            circuit[c] = 1 # the circuit originally does not contain 
            interface.append(circuit)
            
        assert self.connectivity(separation) == len(interface)
        return interface
    
    def similarityMatrix(self, stype):
        """ computes a similarity matrix on the elements
        
        For each pair of elements the similarity matrix has an entry that
        gives the similarity.
        
        stype gives what kind of similarity matrix should be computed.
        Possibilities are:
        direct - direct from the columns of the original coefficient matrix 
            (needs the original coefficient matrix and in the case of 
            contractions the result is not well defined)
        __cols - normalized from the columns (rows are orthogonal)
        rows - normalized from the rows (columns are orthogonal)
        combined - takes max(__cols, rows)
        """
        if stype=="direct":
            return self.__directSimilarityMatrix()
        elif stype=="cols":
            return self.__colSimilarityMatrix()
        elif stype=="rows":
            return self.__rowSimilarityMatrix()
        elif stype=="combined":
            rowSim = self.__rowSimilarityMatrix()
            colSim = self.__colSimilarityMatrix()
            return rowSim.max(colSim)
        else:
            print('unknown option')
            assert False
        
    def __directSimilarityMatrix(self):
        mat = self.matrix.toNumpy()
        
        elems = self.elems
                        
        # compute similarity matrix
        similarity = numpy.zeros((len(elems),len(elems)))
        # precompute norms
        norms = numpy.zeros(len(elems))
        for i in range(len(elems)):
            norms[i] = numpy.sqrt(numpy.dot(mat[:,i],mat[:,i]))
        printMod = ceil(1.0e9 / (len(elems)**2))
        for i in range(len(elems)):
            if i % printMod == 0:
                print('building similarity matrix: %d of %d' % (i, len(elems)))
            for j in range(len(elems)):
                if i != j:
                    similarity[i,j] = numpy.abs(numpy.dot(mat[:,j],mat[:,i]))
                    similarity[i,j] /= norms[i] * norms[j]
                else:
                    similarity[i,j] = -1
                
        return SimilarityMatrix(similarity, elems)
        
    def __rowSimilarityMatrix(self):
        mat = self.rep.toNumpy()
        #mat = sympy.matrix2numpy(self.rep)
        mat = numpy.concatenate((numpy.identity(len(self.nbasis)), mat))
        
        elems = self.nbasis + self.basis
        
        printMod = ceil(1.0e9 / (len(elems)**2))
        # GramSchmidt orthonornalization
        for i in range(len(self.nbasis)):
            if i % printMod == 0 and i > 0:
                print('normalizing matrix: %d of %d' % (i, len(elems)))
            for j in range(i-1):
                mat[:,i] -= numpy.dot(mat[:,j],mat[:,i]) * mat[:,j]
                assert numpy.abs(numpy.dot(mat[:,i],mat[:,j])) < 1e-7
            mat[:,i] /= numpy.sqrt(numpy.dot(mat[:,i],mat[:,i]))
            assert abs(numpy.dot(mat[:,i],mat[:,i]) - 1) < 1e-7
                
        # compute similarity matrix
        similarity = numpy.zeros((len(elems),len(elems)))
        # precompute norms
        norms = numpy.zeros(len(elems))
        for i in range(len(elems)):
            norms[i] = numpy.sqrt(numpy.dot(mat[i,:],mat[i,:]))
        for i in range(len(elems)):
            if i % printMod == 0 and i > 0:
                print('building similarity matrix: %d of %d' % (i, len(elems)))
            for j in range(len(elems)):
                if i != j:
                    similarity[i,j] = numpy.abs(numpy.dot(mat[j,:],mat[i,:]))
                    similarity[i,j] /= norms[i] * norms[j]
                else:
                    similarity[i,j] = -1
                
        return SimilarityMatrix(similarity, elems)
        
    def __colSimilarityMatrix(self):
        # we cannot just reuse __rowSimilarityMatrix, because that would
        # give a different order of the elements
        mat = self.rep.toNumpy()
        #mat = sympy.matrix2numpy(self.rep)
        mat = numpy.concatenate((mat, numpy.identity(len(self.basis))),1)
        
        elems = self.nbasis + self.basis
        
        printMod = ceil(1.0e9 / (len(elems)**2))
        # GramSchmidt orthonornalization
        for i in range(len(self.basis)):
            if i % printMod == 0 and i > 0:
                print('normalizing matrix: %d of %d' % (i, len(elems)))
            for j in range(i-1):
                mat[i,:] -= numpy.dot(mat[j,:],mat[i,:]) * mat[j,:]
            mat[i,:] /= numpy.sqrt(numpy.dot(mat[i,:],mat[i,:]))
        
        # compute similarity matrix
        similarity = numpy.zeros((len(elems),len(elems)))
        # precompute norms
        norms = numpy.zeros(len(elems))
        for i in range(len(elems)):
            norms[i] = numpy.sqrt(numpy.dot(mat[:,i],mat[:,i]))
        for i in range(len(elems)):
            if i % printMod == 0 and i > 0:
                print('building similarity matrix: %d of %d' % (i, len(elems)))
            for j in range(len(elems)):
                if i != j:
                    similarity[i,j] = numpy.abs(numpy.dot(mat[:,j],mat[:,i]))
                    similarity[i,j] /= norms[i] * norms[j]
                else:
                    similarity[i,j] = -1
                
        return SimilarityMatrix(similarity, elems)
        
    
def fromMatrix(matrix, cols):
    """ builds the reduced representation matrix
    
        Given a matrix of columns (self.matrix) it builds the reduced
        representation matrix (self.rep) that has rows for each basis 
        element and columns for each non-basis element.
        
        From this form the fundamental cycles can easily be read of and
        other operations like deletion and contraction are easy.
    """
    assert matrix.cols() <= len(cols), 'for each column a name is needed'
    inv_list = [] # basis vectors in (transformed) triangular form
    row_perm = [] # stores which index we use for each basis element
    # inv_cycle is matrix to reconstruct basis vectors, 
    # used to compute cycles for non-basis elements
    #inv_cycle = sympy.zeros(len(cols)) 
    inv_cycle = Matrix()
    basis = [] # names of basis elements
    nbasis = [] # names of non-basis elements
    repRawData = [] # raw form of representation matrix
    
    time_build_inv_cycle_col = 0
    time_find_row = 0
    time_add_dependent = 0
    time_add_basis = 0
    
    if len(cols) >= 100:
        print("creating matroid representation: ", end="")
    for i in range(len(cols)):
        if i > 0 and i % 100 == 0:
            print(i, end=" ")
        c = matrix[:,i]
        # eliminate existing coordinates (in row_perm)
        stime = time.time()
        res = c.copy() # residual
        # in case its linearly independent, store the linear transformation
        #inv_cycle_col = sympy.zeros(len(cols),c=1)
        inv_cycle_col = Matrix()
        inv_cycle_col[len(row_perm),0] = 1
        for j in range(len(row_perm)):
            coef = c[row_perm[j],0]
            if coef != 0:
                # update residual
                res.coladd(0,inv_list[j], -coef)
                #res -= inv_list[j]*coef
                inv_cycle_col.coladd(0, inv_cycle[:,j], -coef)
                #inv_cycle_col -= inv_cycle[:,j]*coef
            
        time_build_inv_cycle_col += time.time() - stime
        stime = time.time()    
        # check if res contains non-zero entries and
        # if yes, find entry with minimum no. of non-zeros in other __cols
                
        min_index = -1
        min_count = -1
        for j in range(matrix.rows()):
            if res[j,0] != 0:
                #print(row_perm)
                #print(res)
                assert j not in row_perm
                # count no. non-zeros in other __cols
                count = 0
                for x in inv_list:
                    if x[j,0] != 0:
                        count += 1
                # update min_count and min_index
                if min_index < 0 or count < min_count:
                    min_index = j
                    min_count = count
       
        time_find_row += time.time() - stime
        stime = time.time() 
        if min_index == -1:
            # the column i is linearly dependent
            # therefore, add fundamental circuit and add non-basis element
            
            # reconstruct linear combination of original vectors (cycle)
            sub_inv = inv_cycle[0:len(basis),0:len(basis)]
#            v = sub_inv*c.extract(row_perm,[0])
            v = sub_inv*c[row_perm,[0]]
            
            # new non-basis element
            nbasis.append(cols[i])
            # corresponding cycle
            repRawData.append(-v)
            time_add_dependent += time.time() - stime
        else:
            # the column i is linearly independent of other columns
            # therefore, add to basis
            basis.append(cols[i])
            
            # normalize such that the new basis-vector has entry 1 in row j
            # we choose min_index as the 1-entry such that we minimize the
            # number of transformations
            div = res[min_index,0] 
            inv_cycle_col /= div
            # store transformation
            inv_cycle.coladd(len(basis)-1, inv_cycle_col)
            #inv_cycle[:,len(basis)-1] = inv_cycle_col
            res /= div
            
            # eliminate non-zero coefficients of existing basis-vectors 
            # in row min_index
            for k in range(len(inv_list)):
                x = inv_list[k]
                if x[min_index,0] != 0:
                    # also update the transformation matrix
                    inv_cycle.coladd(k, inv_cycle_col, -x[min_index,0])
                    #inv_cycle[:,k] -= inv_cycle_col*x[min_index,0]
                    x.coladd(0, res, -x[min_index,0])
                    #inv_list[k] -= res*x[min_index,0]
            
            # add basis-vector
            inv_list.append(res)
            row_perm.append(min_index)
            time_add_basis += time.time() - stime
   
    if len(cols) >= 100:
        print("") # to finish the line
        
#    print("time build inv cycle col: %s" % time_build_inv_cycle_col)
#    print("time find row %s" % time_find_row)
#    print("time add dependent %s" % time_add_dependent)
#    print("time add basis %s" % time_add_basis)
    # finally, build representation matrix
    #rep = sympy.zeros(len(basis), len(nbasis))
    rep = Matrix()
    for i in range(len(repRawData)):
        rep[:,i] = repRawData[i] # TODO: I think this can be greatly improved!
       
    m = Matroid(rep, basis, nbasis)
    m.matrix = matrix
    m.elems = cols # just make sure that the order is correct
    #clear_cache()
    return m

class SimilarityMatrix:
    
    def __init__(self, matrix, elems):
        """ builds a similarity matrix out of a matrix and element names.
        
        the matrix must be square and have as many rows/columns as elems has
        entries
        """
        self.matrix = matrix
        self.elems = elems
    
    def max(self, simmatrix):
        assert simmatrix.elems == self.elems
        mat = numpy.maximum(self.matrix, simmatrix.matrix)
        
        return SimilarityMatrix(mat, self.elems)
    
    def __getitem__(self, i):
        return self.matrix[i]
    
    def __setitem__(self, i, v):
        self.matrix[i] = v