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

TEMPDATA = {}

class TestReagents(object):
    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.mcore = cbmpy.readSBML3FBC('cbmpy_test_core')
        klass.sboterm1 = 'SBO:1234567'
        if 'HOME' in os.environ:
            klass.CDIR = os.environ['HOME']
        else:
            klass.CDIR = os.getcwd()

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        del klass.mcore

    def roundTripModelV1(self, model):
        modname = 'rtv1.xml'
        cbmpy.writeSBML3FBC(model, modname, directory=self.CDIR)
        del model
        model = cbmpy.readSBML3FBC(modname, work_dir=self.CDIR)
        return model

    def test_param_create(self):
        M = self.mcore.clone()

        pid = 'par1'
        pval = 1.0
        pname = 'parameter1'
        par1 = cbmpy.CBModel.Parameter(pid, pval, pname, True)
        M.addParameter(par1)
        pout = M.getParameter(pid)
        assert pout.getId() == pid
        assert pout.getName() == pname
        assert pout.getValue() == pval
        assert pout.constant == True
        del M

    def test_set_coefficient(self):
        M = self.mcore.clone()

        pval = 10.0
        rid = 'R01'
        suid = 'X0'
        prid = 'A'
        R = M.getReaction(rid)
        sub = R.getReagentWithSpeciesRef(suid)
        prod = R.getReagentWithSpeciesRef(prid)

        sub.setCoefficient(-pval)
        prod.setCoefficient(pval)

        #M = self.roundTripModelV1(M)

        R = M.getReaction(rid)
        sub = R.getReagentWithSpeciesRef(suid)
        prod = R.getReagentWithSpeciesRef(prid)

        assert sub.getCoefficient() == -pval
        assert prod.getCoefficient() == pval

        del M

    def test_set_par_coefficient(self):
        M = self.mcore.clone()

        pname = 'parameter1'
        pval = 10.0
        rid = 'R01'
        suid = 'X0'
        prid = 'A'

        par_sub = cbmpy.CBModel.Parameter('par1', -pval, pname, True)
        par_prod = cbmpy.CBModel.Parameter('par2', pval, pname, True)
        M.addParameter(par_sub)
        M.addParameter(par_prod)

        R = M.getReaction(rid)
        sub = R.getReagentWithSpeciesRef(suid)
        prod = R.getReagentWithSpeciesRef(prid)

        sub.setCoefficient(par_sub)
        prod.setCoefficient(par_prod)

        #M = self.roundTripModelV1(M)

        R = M.getReaction(rid)
        sub = R.getReagentWithSpeciesRef(suid)
        prod = R.getReagentWithSpeciesRef(prid)

        assert sub.getCoefficient() == -pval
        assert prod.getCoefficient() == pval

        del M







    #def test_set_coefficient(self):
        #M = self.mcore.clone()
        #for c in M.reactions:
            #c.setSBOterm(self.sboterm1)
        #M = self.roundTripModelV1(M)
        #for c in M.reactions:
            #assert self.sboterm1 == c.getSBOterm()
        #del M
