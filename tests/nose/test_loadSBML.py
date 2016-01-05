import os, sys, time, zipfile, json
cDir = os.path.dirname(os.path.abspath(__file__))
MDIR = os.path.join(os.path.split(cDir)[0], 'models')
MDIR = cDir

import nose
from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises

import cbmpy

F = open(os.path.join(MDIR, 'data.json'))
DATA = json.load(F)
F.close()

class TestReadSBML(object):
    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_load_L3FBCV2_iJR904(self):
        cmod = cbmpy.readSBML3FBC(os.path.join(MDIR, 'L3FBCV2_iJR904.glc.xml'))
        assert_not_equal(cmod, None)
        
    def test_load_L3FBCV1_iJR904(self):
        cmod = cbmpy.readSBML3FBC(os.path.join(MDIR, 'L3FBCV1_iJR904.glc.xml'))
        assert_not_equal(cmod, None)

    def test_load_L2FBA_iJR904(self):
        cmod = cbmpy.readSBML2FBA(os.path.join(MDIR, 'L2FBA_iJR904.glc.xml'))
        assert_not_equal(cmod, None)

    def test_load_L2CBR_iJR904(self):
        cmod = cbmpy.readCOBRASBML(os.path.join(MDIR, 'L2CBR_iJR904.glc.xml'),\
                                                delete_intermediate=True)
        assert_not_equal(cmod, None)


class TestRunFBA(object):
    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.m = {}
        for m in DATA:
            if m.startswith('L2CBR_'):
                cmod = cbmpy.readCOBRASBML(os.path.join(MDIR, m),\
                                           delete_intermediate=True)
            elif m.startswith('L2FBA_'):
                cmod = cbmpy.readSBML2FBA(os.path.join(MDIR, m))
            elif m.startswith('L3FBCV1_') or m.startswith('L3FBCV2_'):
                cmod = cbmpy.readSBML3FBC(os.path.join(MDIR, m))
            klass.m[m] = cmod

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        del klass.m

    def test_load(self):
        for m in self.m:
            assert_not_equal(self.m[m], None)
            
    def test_run_fba_cobra(self):
        cmod = cbmpy.readCOBRASBML(os.path.join(MDIR, 'L2CBR_iJR904.glc.xml'),\
                                                        delete_intermediate=True)
        cmod2 = cbmpy.readSBML3FBC(os.path.join(MDIR, 'L3FBCV1_iJR904.glc.xml'))
        fba1 = cbmpy.doFBA(cmod)
        fba2 = cbmpy.doFBA(cmod2)
        del cmod, cmod2
        assert_almost_equal(fba1, fba2)
        
    def test_run_L3FBCV2_iJR904(self):
        m = 'L3FBCV2_iJR904.glc.xml'
        cmod = self.m[m].clone()
        res = cbmpy.doFBA(cmod)
        del cmod
        assert_almost_equal(DATA[m]['FBA'], res)
        
    def test_run_L3FBCV1_iJR904(self):
        m = 'L3FBCV1_iJR904.glc.xml'
        cmod = self.m[m].clone()
        res = cbmpy.doFBA(cmod)
        del cmod
        assert_almost_equal(DATA[m]['FBA'], res)
        
    def test_run_L2FBA_iJR904(self):
        m = 'L2FBA_iJR904.glc.xml'
        cmod = self.m[m].clone()
        res = cbmpy.doFBA(cmod)
        del cmod
        assert_almost_equal(DATA[m]['FBA'], res)


class TestBiGG2(object):
    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.m = {}
        for m in DATA:
            if m.startswith('BIGG2_'):
                cmod = cbmpy.readSBML3FBC(os.path.join(MDIR, m))
                klass.m[m] = cmod

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        del klass.m

    def test_load(self):
        for m in self.m:
            assert_not_equal(self.m[m], None)
            
    def test_run_BIGG2_iIT341(self):
        m = 'BIGG2_iIT341.xml'
        cmod = self.m[m].clone()
        res = cbmpy.doFBA(cmod)
        del cmod
        assert_almost_equal(DATA[m]['FBA'], res)
        
    def test_run_BIGG2_iAF1260(self):
        m = 'BIGG2_iAF1260.xml'
        cmod = self.m[m].clone()
        res = cbmpy.doFBA(cmod)
        del cmod
        assert_almost_equal(DATA[m]['FBA'], res)

        