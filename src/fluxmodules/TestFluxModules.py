"""
CBMPy: fluxmodules main functionality test module

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

import cbmpy as cbm

from . import fluxmodules

class TestToy(unittest.TestCase):
    """ uses a simple toy model to test correctness """

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "toy_model_maarleveld2014.l3.xml"
        
        self.cmod = cbm.CBRead.readSBML3FBC(model_file, model_dir)


    def tearDown(self):
        pass


    def testComputeModules(self):
        cbm.CBSolver.FluxVariabilityAnalysis(self.cmod)  # @UndefinedVariable
        mods = fluxmodules.computeModules(self.cmod)
        self.assertEqual(len(mods),3)
        
        # expected results
        expected = {frozenset(['R03','R02','R04']),
                    frozenset(['R07','R06','R09','R08','R10']),
                    frozenset(['R12','R13','R14','R15'])}
        expectedVariable = set()
        for e in expected:
            expectedVariable.update(e)
        
        # check if every computed module is indeed an expected module
        variable = set()
        for m in mods:
            self.assertTrue(frozenset(m) in expected)
            variable.update(m)
        pass
        
        # to make sure that we did not miss any module, check if the set
        # of variable reactions is the same
        self.assertSetEqual(variable, expectedVariable)

class TestEcoli(unittest.TestCase):
    """ uses a genome-scale metabolic networks """

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "Ecoli_iAF1260_ox.glc.cobra.xml"
        
        self.cmod = cbm.CBRead.readCOBRASBML(model_file, model_dir)


    def tearDown(self):
        pass


    def testComputeModules(self):
        cbm.CBSolver.FluxVariabilityAnalysis(self.cmod)  # @UndefinedVariable
        mods = fluxmodules.computeModules(self.cmod)
        self.assertEqual(len(mods),13)
        
        # expected resultsr_
        expected = {frozenset(['R_GLUt4pp', 'R_CAt6pp', 'R_ACt2rpp', 'R_ACt4pp', 
                               'R_GLUABUTt7pp', 'R_GLYCLTt2rpp', 'R_PROt2rpp', 
                               'R_NAt3pp', 'R_ABUTt2pp', 'R_CA2t3pp', 
                               'R_SERt2rpp', 'R_SERt4pp', 'R_GLUt2rpp', 
                               'R_GLYCLTt4pp', 'R_THRt4pp', 'R_THRt2rpp', 
                               'R_PROt4pp']),
                    frozenset(['R_EX_h_e_', 'R_FE2tex', 'R_Htex', 'R_EX_h2o_e_',
                               'R_FE3tex', 'R_EX_o2_e_', 'R_EX_fe3_e_', 
                               'R_EX_fe2_e_', 'R_O2tex', 'R_FEROpp', 
                               'R_H2Otex']),
                    frozenset(['R_PRPPS', 'R_RNDR4', 'R_RNDR2', 'R_RNDR3', 
                               'R_GLCP', 'R_PPKr', 'R_GRXR', 'R_R15BPK', 
                               'R_PUNP1', 'R_RNDR4b', 'R_RNTR3c', 'R_PAPSR', 
                               'R_FLDR', 'R_ADK3', 'R_ACKr', 'R_ADNK1', 
                               'R_RNTR2c', 'R_GLDBRAN2', 'R_PPM', 'R_RNDR2b', 
                               'R_RNDR1b', 'R_NDPK7', 'R_NDPK6', 'R_NDPK5', 
                               'R_NDPK3', 'R_NDPK2', 'R_GLBRAN2', 'R_R1PK', 
                               'R_NDPK8', 'R_GLCP2', 'R_PTAr', 'R_NDPK1', 
                               'R_GLGC', 'R_RNTR4c', 'R_RNDR1', 'R_GLCS1', 
                               'R_GTHOr', 'R_TRDR', 'R_RNDR3b', 'R_ADK1', 
                               'R_ACS', 'R_PAPSR2', 'R_ADPT', 'R_RNTR1c']),
                    frozenset(['R_DMPPS', 'R_IPDPS', 'R_IPDDI']),
                    frozenset(['R_ACCOAL', 'R_SUCOAS', 'R_PPAKr', 
                               'R_PPCSCT', 'R_PTA2']),
                    frozenset(['R_CRNDt2rpp', 'R_CRNt8pp', 'R_CRNt2rpp']),
                    frozenset(['R_GLYCTO2', 'R_GLYCTO3', 'R_GLYCTO4', 'R_ASPO3',
                               'R_ASPO5', 'R_ASPO4', 'R_DHORD5', 'R_SUCDi',
                               'R_DHORD2', 'R_FRD3', 'R_FRD2', 'R_NADH17pp', 
                               'R_NADH18pp', 'R_NADH16pp']),
                    frozenset(['R_TRSARr', 'R_HPYRI', 'R_HPYRRx']),
                    frozenset(['R_ACACT1r', 'R_KAT1']),
                    frozenset(['R_PFK', 'R_PYK', 'R_FBA', 'R_DHAPT', 'R_F6PA']),
                    frozenset(['R_KAS15', 'R_ACOATA', 'R_KAS14']),
                    frozenset(['R_INDOLEt2rpp', 'R_INDOLEt2pp']),
                    frozenset(['R_VPAMT', 'R_VALTA', 'R_ALATA_L'])
                   }
        expectedVariable = set()
        for e in expected:
            expectedVariable.update(e)
        
        # check if every computed module is indeed an expected module
        variable = set()
        for m in mods:
            self.assertTrue(frozenset(m) in expected)
            variable.update(m)
        pass
        
        # to make sure that we did not miss any module, check if the set
        # of variable reactions is the same
        self.assertSetEqual(variable, expectedVariable)

class TestITM686(unittest.TestCase):
    """ uses a simple toy model to test correctness """

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "iTM686_sbml3.xml"
        
        self.cmod = cbm.CBRead.readSBML3FBC(model_file, model_dir)


    def tearDown(self):
        pass


    def testComputeModules(self):
        cbm.CBSolver.FluxVariabilityAnalysis(self.cmod)  # @UndefinedVariable
        mods = fluxmodules.computeModules(self.cmod)
        print(mods)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()