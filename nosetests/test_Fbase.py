import os, sys, time, zipfile, json
cDir = os.path.dirname(os.path.abspath(__file__))
import cbmpy
MDIR = cDir

import nose
from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises



F = open(os.path.join(MDIR, 'data.json'))
DATA = json.load(F)
F.close()

TEMPDATA = {}

class TestBase:
    def test_import(self):
        import cbmpy

class TestSBO(object):
    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.mcore = cbmpy.readSBML3FBC('cbmpy_test_core')
        klass.sboterm1 = 'SBO:1234567'
        klass.cDir = cDir

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        del klass.mcore

    def roundTripModelV1(self, model):
        modname = 'rtv1.xml'
        cbmpy.writeSBML3FBCV2(model, modname, directory=self.cDir)
        del model
        model = cbmpy.readSBML3FBC(os.path.join(self.cDir, modname), xoptions={'validate' : True})
        return model

    def test_sbo_reactions_getset(self):
        M = self.mcore.clone()
        for c in M.reactions:
            c.setSBOterm(self.sboterm1)
        for c in M.reactions:
            print(self.sboterm1, c.getSBOterm())
            assert self.sboterm1 == c.getSBOterm()
        del M

    def test_sbo_species_getset(self):
        M = self.mcore.clone()
        for c in M.species:
            c.setSBOterm(self.sboterm1)
        for c in M.species:
            assert self.sboterm1 == c.getSBOterm()
        del M

    def test_sbo_compartments_getset(self):
        M = self.mcore.clone()
        for c in M.compartments:
            c.setSBOterm(self.sboterm1)
        for c in M.compartments:
            assert self.sboterm1 == c.getSBOterm()
        del M

    def test_sbo_reactions_getset_wrt(self):
        M = self.mcore.clone()
        for c in M.reactions:
            c.setSBOterm(self.sboterm1)
        M = self.roundTripModelV1(M)
        for c in M.reactions:
            assert self.sboterm1 == c.getSBOterm()
        del M

    def test_sbo_species_getset_wrt(self):
        M = self.mcore.clone()
        for c in M.species:
            c.setSBOterm(self.sboterm1)
        M = self.roundTripModelV1(M)
        for c in M.species:
            assert self.sboterm1 == c.getSBOterm()
        del M

    def test_sbo_compartments_getset_wrt(self):
        M = self.mcore.clone()
        for c in M.compartments:
            c.setSBOterm(self.sboterm1)
        M = self.roundTripModelV1(M)
        for c in M.compartments:
            assert self.sboterm1 == c.getSBOterm()
        del M
