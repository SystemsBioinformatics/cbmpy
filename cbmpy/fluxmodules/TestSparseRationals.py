"""
CBMPy: fluxmodules sparse rationals test module
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

import unittest

from sparserationals import Matrix

from fractions import Fraction

class Test(unittest.TestCase):


    def setUp(self):
        self.mat = Matrix([[-2, -1,  1,  0,  0,  0,  0,  0,  0],
                           [ 1,  0,  0, -1,  0,  0,  0,  0,  0],
                           [ 0,  1, -1,  1, -1, -1,  0,  0,  0],
                           [ 1,  0,  0,  0,  1,  0, -1,  0,  0],
                           [ 0,  0,  0,  0,  0,  1,  1, -1,  1],
                           [ 0,  0,  0,  0,  0,  0,  0,  1, -1]])


    def tearDown(self):
        pass


    def testPrintNumpy(self):
        self.assertEquals(self.mat.cols(), 9)
        self.assertEquals(self.mat.rows(), 6)
        print(self.mat.toNumpy())

    def testValueAccess(self):
        self.assertEqual(self.mat[0,0], -2)
        self.assertEqual(self.mat[5,3], 0)
        self.assertEquals(self.mat[3,6], -1)
        
    def testSubMatrix(self):
        submat = self.mat[:,[2,4]]
        self.assertEqual(submat[0,0], 1)
        self.assertEqual(submat[1,0], 0)
        self.assertEqual(submat[2,0], -1)
        self.assertEqual(submat[3,0], 0)
        self.assertEqual(submat[4,0], 0)
        self.assertEqual(submat[5,0], 0)
        self.assertEqual(submat[0,1], 0)
        self.assertEqual(submat[1,1], 0)
        self.assertEqual(submat[2,1], -1)
        self.assertEqual(submat[3,1], 1)
        self.assertEqual(submat[4,1], 0)
        self.assertEqual(submat[5,1], 0)
        
    def testTranspose(self):
        t = self.mat.transpose()
        tt = t.transpose()
        tt -= self.mat
        for i in range(tt.rows()):
            for j in range(tt.cols()):
                self.assertEqual(tt[i,j], 0)
    
    def testSetItem(self):
        submat = self.mat[[1,3],1:5]
        c = self.mat.copy()
        submat *= Fraction(-5)
        c[[1,3],1:5] += submat / 5 
        
        for i in [1,3]:
            for j in range(*(1,5)):
                self.assertEqual(c[i,j], 0)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()