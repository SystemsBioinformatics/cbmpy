"""
CBMPy: fluxmodules decomposiiton test module

This module only exists for testing purposes
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

import unittest
import os

from pyscescbm import CBRead

from . import matroid
from .decomposition import Decomposition, Vertex
from .sparserationals import Matrix

class TestDecomposition(unittest.TestCase):
 
 
    def setUp(self):
        """ create test matroid """
        mat = Matrix([[-2, -1,  1,  0,  0,  0,  0,  0,  0],
                      [ 1,  0,  0, -1,  0,  0,  0,  0,  0],
                      [ 0,  1, -1,  1, -1, -1,  0,  0,  0],
                      [ 1,  0,  0,  0,  1,  0, -1,  0,  0],
                      [ 0,  0,  0,  0,  0,  1,  1, -1,  1],
                      [ 0,  0,  0,  0,  0,  0,  0,  1, -1]])
         
        elems = [1,2,3,4,5,6,7,8,9]
        self.matroid = matroid.fromMatrix(mat, elems)
 
        self.decomp = Decomposition(self.matroid)
 
 
    def tearDown(self):
        pass
 
 
    def testZeroModules(self):
        self.decomp.splitZeroModules(self.decomp.root)
         
        self.assertEqual(len(self.decomp.root.edges),2)
        for e in self.decomp.root.edges.values():
            self.assertEqual(len(e), 0)
         
        for e in self.decomp.root.edges.keys():
            self.assertIsInstance(e, Vertex)
            self.assertEqual(len(e.edges), len(e.matroid.elems)+1)
            if len(e.matroid.elems) == 2:
                self.assertSetEqual(set(e.matroid.elems), {8,9})
            else:
                self.assertSetEqual(set(e.matroid.elems), {1,2,3,4,5,6,7})
 
    def testOneModules(self):
        for v in self.decomp.listNonLeaves():
            self.decomp.splitOneModules(v)
         
        for v in self.decomp.listNonLeaves():
            print(v)
            print(v.edges)
            print("")
 
        self.assertTrue(self.decomp.isParallelSeries())
         
        for v in self.decomp.listNonLeaves():
            for (w, interface) in v.edges.items():
                sep = w.getLeaves(v)
                print(sep)
                conn = self.matroid.connectivity(sep)
                self.assertLessEqual(conn, 1)
                self.assertEqual(len(interface), conn)
                 
class TestDecomposition2(unittest.TestCase):
    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "toy_model_maarleveld2014.l3.xml"
          
        self.cmod = CBRead.readSBML3FBC(model_file, model_dir)
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
  
        print(self.matroid.basis)
        print(self.matroid.nbasis)
        print(self.matroid.rep.toNumpy())
  
    def tearDown(self):
        pass
  
    def testPoolmanDirectRows(self):
        decomp = Decomposition(self.matroid)
        sim = self.matroid.similarityMatrix("rows")
#        print("rows:")
#        print(sim.matrix)
        # we expect coparallel elements to have similarity 1
        coparallel = self.matroid.coparallel()
        for group in coparallel:
            for e in group:
                for f in group:
                    if e != f:
                        ei = sim.elems.index(e)
                        fi = sim.elems.index(f)
                        self.assertAlmostEqual(sim[ei,fi], 1, 
                                msg="at %s(%d) %s(%d)" % (e, ei, f, fi))
                    
        decomp.poolmanMethod(decomp.root, sim)
#        for v in decomp.listNonLeaves():
#            print(v)
#            print(v.edges)
#            print("")
        write(decomp.makeGraphViz(),"poolman_direct_rows.gv")
  
    def testPoolmanDirectCols(self):
        decomp = Decomposition(self.matroid)
        sim = self.matroid.similarityMatrix("cols")
#        print("cols:")
#        print(sim.matrix)
        # we expect parallel elements to have similarity 1
        parallel = self.matroid.parallel()
        for group in parallel:
            for e in group:
                for f in group:
                    if e != f:
                        ei = sim.elems.index(e)
                        fi = sim.elems.index(f)
                        self.assertAlmostEqual(sim[ei,fi], 1, 
                                msg="at %s(%d) %s(%d)" % (e, ei, f, fi))
                        
        decomp.poolmanMethod(decomp.root, sim)
#        for v in decomp.listNonLeaves():
#            print(v)
#            print(v.edges)
#            print("")
        write(decomp.makeGraphViz(),"poolman_direct_cols.gv")
 
    def testPoolmanDirectDirect(self):
        decomp = Decomposition(self.matroid)
        sim = self.matroid.similarityMatrix("direct")
#        print("cols:")
#        print(sim.matrix)
        # we expect parallel elements to have similarity 1
        parallel = self.matroid.parallel()
        for group in parallel:
            for e in group:
                for f in group:
                    if e != f:
                        ei = sim.elems.index(e)
                        fi = sim.elems.index(f)
                        self.assertAlmostEqual(sim[ei,fi], 1, 
                                msg="at %s(%d) %s(%d)" % (e, ei, f, fi))
                        
        decomp.poolmanMethod(decomp.root, sim)
#        for v in decomp.listNonLeaves():
#            print(v)
#            print(v.edges)
#            print("")
        write(decomp.makeGraphViz(),"poolman_direct_direct.gv")
 
    def testPoolmanDirectCombinded(self):
        decomp = Decomposition(self.matroid)
        sim = self.matroid.similarityMatrix("combined")
#        print("combined:")
#        print(sim.matrix)
        # we expect coparallel elements to have similarity 1
        coparallel = self.matroid.coparallel()
        for group in coparallel:
            for e in group:
                for f in group:
                    if e != f:
                        ei = sim.elems.index(e)
                        fi = sim.elems.index(f)
                        self.assertAlmostEqual(sim[ei,fi], 1, 
                                msg="at %s(%d) %s(%d)" % (e, ei, f, fi))
        # we expect parallel elements to have similarity 1
        parallel = self.matroid.parallel()
        for group in parallel:
            for e in group:
                for f in group:
                    if e != f:
                        ei = sim.elems.index(e)
                        fi = sim.elems.index(f)
                        self.assertAlmostEqual(sim[ei,fi], 1, 
                                msg="at %s(%d) %s(%d)" % (e, ei, f, fi))
                        
        decomp.poolmanMethod(decomp.root, sim)
#        for v in decomp.listNonLeaves():
#            print(v)
#            print(v.edges)
#            print("")
        write(decomp.makeGraphViz(),"poolman_direct_combined.gv")

    def testOneModules(self):
        decomp = Decomposition(self.matroid)
        decomp.splitOneModules(decomp.root)
        write(decomp.makeGraphViz(),"one_modules.gv")

    def testPoolmanRows(self):
        decomp = Decomposition(self.matroid)
        decomp.splitOneModules(decomp.root)
        sim = decomp.root.matroid.similarityMatrix("rows")
        decomp.poolmanMethod(decomp.root, sim)
        write(decomp.makeGraphViz(),"poolman_rows.gv")


#@unittest.skip("takes quite long")
class TestDecompositionITM686(unittest.TestCase):
    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "iTM686_sbml3.xml"
          
        self.cmod = CBRead.readSBML3FBC(model_file, model_dir)
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
  
    def tearDown(self):
        pass
    
    @unittest.skip("takes quite long")    
    def testPoolmanRows(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("rows")
                decomp.poolmanMethod(v, sim)
           
        write(decomp.makeGraphViz(),"iTM686_poolman_rows.gv")
        self.assertTrue(decomp.verifyFullyBranched())
        self.assertTrue(decomp.verifyEdgeWidth())
        print("width rows = %d" % decomp.getWidth())
        
    @unittest.skip("takes quite long")    
    def testPoolmanCols(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("cols")
                decomp.poolmanMethod(v, sim)
           
        write(decomp.makeGraphViz(),"iTM686_poolman_cols.gv")
        #self.assertTrue(decomp.verifyFullyBranched())
        #self.assertTrue(decomp.verifyEdgeWidth())
        print("width cols = %d" % decomp.getWidth())
    
    #@unittest.skip("takes quite long")    
    def testPoolmanCombined(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("combined")
                decomp.poolmanMethod(v, sim)
          
        write(decomp.makeGraphViz(),"iTM686_poolman_combined.gv")
        #self.assertTrue(decomp.verifyFullyBranched())
        #self.assertTrue(decomp.verifyEdgeWidth())
        print("width combined = %d" % decomp.getWidth())
        
    @unittest.skip("takes quite long")    
    def testPoolmanDirect(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        sim = decomp.matroid.similarityMatrix("direct")
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                decomp.poolmanMethod(v, sim)
          
        write(decomp.makeGraphViz(),"iTM686_poolman_direct.gv")
        #self.assertTrue(decomp.verifyFullyBranched())
        #self.assertTrue(decomp.verifyEdgeWidth())
        print("width direct = %d" % decomp.getWidth())
        
        
class TestDecompositionEcore(unittest.TestCase):
    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "ecoli_core_COBRA.xml"
          
        self.cmod = CBRead.readCOBRASBML(model_file, model_dir)
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
  
    def tearDown(self):
        pass
    
    def testPoolmanRows(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("rows")
                decomp.poolmanMethod(v, sim)
         
        write(decomp.makeGraphViz(cmod=self.cmod),"ecoli_core_poolman_rows.gv")
        self.assertTrue(decomp.verifyFullyBranched())
        self.assertTrue(decomp.verifyEdgeWidth())
        print("width rows = %d" % decomp.getWidth())
        
    
    def testPoolmanCols(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("cols")
                decomp.poolmanMethod(v, sim)
         
        write(decomp.makeGraphViz(cmod=self.cmod),"ecoli_core_poolman_cols.gv")
        self.assertTrue(decomp.verifyFullyBranched())
        self.assertTrue(decomp.verifyEdgeWidth())
        print("width cols = %d" % decomp.getWidth())
    
    def testPoolmanCombined(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("combined")
                decomp.poolmanMethod(v, sim)
         
        write(decomp.makeGraphViz(cmod=self.cmod),"ecoli_core_poolman_combined.gv")
        self.assertTrue(decomp.verifyFullyBranched())
        self.assertTrue(decomp.verifyEdgeWidth())
        print("width combined = %d" % decomp.getWidth())
    
    def testPoolmanDirect(self):
        decomp = Decomposition(self.matroid)
        decomp.splitSimple()
        sim = decomp.matroid.similarityMatrix("direct")
        for v in decomp.listNonLeaves():
            if not v.isSimple():
                decomp.poolmanMethod(v, sim)
         
        write(decomp.makeGraphViz(cmod=self.cmod),"ecoli_core_poolman_direct.gv")
        self.assertTrue(decomp.verifyFullyBranched())
        self.assertTrue(decomp.verifyEdgeWidth())
        print("width direct = %d" % decomp.getWidth())
    
    
def write(result, filename):
    cDir = os.path.dirname(__file__)
    result_dir = os.path.join(cDir,'../../results')
    result_file = os.path.join(result_dir, filename)
    fo = open(result_file, mode='w')
    fo.write(result)
    fo.close()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()