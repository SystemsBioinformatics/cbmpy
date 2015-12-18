"""
CBMPy: fluxmodules matroid test module

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

#from sympy import Matrix
from .sparserationals import Matrix

from . import matroid

class TestMatroid(unittest.TestCase):
    """ Test methods in matroid.py """

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

    def tearDown(self):
        """ nothing to do"""
        pass


    def testElems(self):
        self.assertEqual(len(self.matroid.elems), 9)
        self.assertEqual(len(self.matroid.basis)+len(self.matroid.nbasis), 
                          len(self.matroid.elems))
        # every basis must have same size of 5
        self.assertEqual(len(self.matroid.basis), 5)
        
    def testBasis(self):
        print(self.matroid.rep.toNumpy())
        self.assertTrue(self.matroid.isIndependent({1,2,4,6,8}))
        self.assertSetEqual(set(self.matroid.nbasis), {3,5,7,9})
        self.assertSetEqual(set(self.matroid.basis), {1,2,4,6,8})
        
        self.assertTrue(self.matroid.isIndependent({1,3,5,9,7}))
        self.assertTrue(self.matroid.isCoIndependent({3,7,9}))
        
        self.assertFalse(self.matroid.isIndependent({2,3}))
        self.assertTrue(self.matroid.isIndependent({1,4,5}))
        self.assertFalse(self.matroid.isIndependent({8,9}))
        self.assertFalse(self.matroid.isIndependent({1,3,6,7,4}))
    
    def testFundamentalCircuit(self):
        # first set basis
        self.assertTrue(self.matroid.isIndependent({1,2,4,6,8}))
        self.assertSetEqual(self.matroid.fundamentalCircuit(3), {2})
        self.assertSetEqual(self.matroid.fundamentalCircuit(5), {1,2,4})
        self.assertSetEqual(self.matroid.fundamentalCircuit(7), {1,2,4,6})
        self.assertDictEqual(self.matroid.fundamentalCircuitValues(3), {2:1})
        self.assertDictEqual(self.matroid.fundamentalCircuitValues(5), {1:-1, 2:2, 4:-1})
        self.assertDictEqual(self.matroid.fundamentalCircuitValues(7), {1:1, 2:-2, 4:1, 6:-1})

        # if we give basis elements, we get cofundamental circuits
        self.assertSetEqual(self.matroid.fundamentalCircuit(8), {9})
        self.assertSetEqual(self.matroid.fundamentalCircuit(4), {5,7})
 
    def testModules(self):
        modules = self.matroid.findModules()
        self.assertEqual(len(modules), 2)
        for m in modules:
            if len(m.elems) == 2:
                self.assertSetEqual(set(m.elems), {8,9})
            else:
                self.assertSetEqual(set(m.elems), {1,2,3,4,5,6,7})
    
    def testParallel(self):
        parallel = self.matroid.parallel()
        self.assertSetEqual(parallel, {frozenset({2,3}),
                                       frozenset({8,9})})
        
    def testCoParallel(self):
        coparallel = self.matroid.coparallel()
        self.assertSetEqual(coparallel, {frozenset({1,4}),
                                         frozenset({8,9}),
                                         frozenset({6,7})})

    def testDelete(self):
        m = self.matroid.delete({7})
        self.assertFalse(7 in m.elems)
        self.assertTrue(m.isIndependent({6}))
        self.assertFalse(m.isCoIndependent({6}))
        self.assertTrue(m.isIndependent({1,2,4,8}))
        
        self.assertTrue(self.matroid.isCoIndependent({6}))
        
    def testContraction(self):
        m = self.matroid.contract({2,4})
        self.assertFalse(m.isIndependent({3}))
        self.assertFalse(m.isIndependent({1,5}))
        self.assertTrue(m.isIndependent({1}))
        self.assertTrue(m.isIndependent({5}))
        
        parallel = m.contract({6}).parallel()
        self.assertSetEqual(parallel, {frozenset({1,5,7}),
                                       frozenset({8,9})})

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()