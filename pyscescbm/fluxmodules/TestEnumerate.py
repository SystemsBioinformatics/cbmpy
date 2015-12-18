'''
Created on Nov 17, 2014

@author: arne
'''
import unittest

import os
from math import isnan

from pyscescbm import CBRead
from pyscescbm import CBTools
from pyscescbm import CBModel

from pyscescbm.solver import CBSolver

from .enumerate import EFMEnumerator
from .sparserationals import Matrix
from .decomposition import Decomposition, Vertex
from . import matroid

class TestEnumerateToy(unittest.TestCase):

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "toy_model_maarleveld2014.l3.xml"
          
        self.cmod = CBRead.readSBML3FBC(model_file, model_dir)
        self.cmod = CBTools.splitReversibleReactions(self.cmod)
        
        r = self.cmod.getReaction('R01')
        print(r.getUpperBound())
        print(r.getLowerBound())
        r.setLowerBound(2.0)
        
        for r in self.cmod.reactions:
            r.setUpperBound(min(1000, r.getUpperBound()))
        
        self.cmod.buildStoichMatrix()
        
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
        
        self.decomp = Decomposition(self.matroid)
        sim = self.matroid.similarityMatrix("combined")
        self.decomp.poolmanMethod(self.decomp.root, sim)


    def tearDown(self):
        pass


    def testEnumerate(self):
        CBSolver.FluxVariabilityAnalysis(self.cmod, objF2constr=False)
        enum = EFMEnumerator(self.cmod, self.decomp.root, None)
        F = enum.enumerateMinimal()
        print(F)
        self.assertEquals(F, set(
            [frozenset(['R01', 'R17', 'R05', 'R10', 'R11', 'R09', 'R18', 'R03_fwd', 'R13_fwd', 'R12', 'R16', 'R14_fwd', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R05', 'R10', 'R11', 'R09', 'R18', 'R02_fwd', 'R12', 'R16', 'R14_fwd', 'R13_fwd']), 
             frozenset(['R01', 'R17', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R03_fwd', 'R13_fwd', 'R12', 'R16', 'R14_fwd', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R02_fwd', 'R12', 'R16', 'R14_fwd', 'R13_fwd']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R10', 'R11', 'R09', 'R18', 'R03_fwd', 'R12', 'R16', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R10', 'R11', 'R09', 'R18', 'R03_fwd', 'R14_rev', 'R13_rev', 'R16', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R10', 'R11', 'R09', 'R18', 'R02_fwd', 'R12', 'R16']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R10', 'R11', 'R09', 'R18', 'R02_fwd', 'R14_rev', 'R16', 'R13_rev']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R03_fwd', 'R12', 'R16', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R03_fwd', 'R14_rev', 'R13_rev', 'R16', 'R04_fwd']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R02_fwd', 'R14_rev', 'R16', 'R13_rev']), 
             frozenset(['R01', 'R17', 'R15', 'R05', 'R07', 'R11', 'R08', 'R06', 'R18', 'R02_fwd', 'R12', 'R16'])]
                                 )
                          )
        evaluateReuse(F, self.decomp.root, None)

@unittest.skip("takes quite long")
class TestEnumerateECore(unittest.TestCase):

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "ecoli_core_COBRA.xml"
          
        self.cmod = CBRead.readCOBRASBML(model_file, model_dir)
        self.cmod = CBTools.splitReversibleReactions(self.cmod)
        for r in self.cmod.reactions:
            r.setUpperBound(min(1000, r.getUpperBound()))
            r.setLowerBound(0)
        
        r = self.cmod.getReaction('R_Biomass_Ecoli_core_w_GAM')
        print(r.getUpperBound())
        print(r.getLowerBound())
        r.setLowerBound(0.1)
        r.setUpperBound(0.1)
        
        self.cmod.buildStoichMatrix()
        
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
        
        self.decomp = Decomposition(self.matroid)
        self.decomp.splitSimple()
        for v in self.decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("combined")
                self.decomp.poolmanMethod(v, sim)

        write(self.decomp.makeGraphViz(cmod=self.cmod),"ecoli_core_split_poolman_combined.gv")
        #self.assertTrue(self.decomp.verifyFullyBranched())
        #self.assertTrue(self.decomp.verifyEdgeWidth())
        print("width = %d" % self.decomp.getWidth())

    def tearDown(self):
        pass


    def testEnumerate(self):
        CBSolver.FluxVariabilityAnalysis(self.cmod, objF2constr=False)
        enum = EFMEnumerator(self.cmod, self.decomp.root, None)
        fObj = CBModel.FluxObjective('R_EX_glc_e_rev_obj', reaction='R_EX_glc_e_rev', coefficient=1)
        objF = CBModel.Objective('minimize_nutrient', 'min')
        objF.addFluxObjective(fObj)
        self.cmod.addObjective(objF, active=True)
        CBSolver.analyzeModel(self.cmod)
        #CBSolver.MinimizeSumOfAbsFluxes(self.cmod, objF2constr=False)
        face = set()
        for r in self.cmod.reactions:
            if r.getValue() > 1e-5:
                face.add(r.getId())
                print('%s: %g' % (r.getId(), r.getValue()) )
        self.assertTrue(enum.test(face))
        F = enum.enumerateMinimal()
        print(F)
        
def evaluateReuse(F, node, parent):
    mod = node.getLeaves(parent)
    Fmod = {}
    for f in F:
        fmod = f.intersection(mod)
        Fmod[fmod] = Fmod.get(fmod, 0) + 1
    FmodVals = Fmod.values()
    FmodVals.sort()
    FmodVals.reverse()
    print(mod)
    print(FmodVals)
    if len(mod) > 10 and isinstance(node, Vertex):
        adj = node.getAdjacent()
        for a in adj:
            if a != parent:
                evaluateReuse(F, a, node)
            

@unittest.skip("takes quite long")
class TestEnumerateEcoliTextbookBigg(unittest.TestCase):

    def setUp(self):
        cDir = os.path.dirname(__file__)
        model_dir = os.path.join(cDir,'../models')
        model_file = "bigg_textbook.xml.l3fbc.xml"
          
        self.cmod = CBRead.readSBML3FBC(model_file, model_dir)
        self.cmod = CBTools.splitReversibleReactions(self.cmod)
        for r in self.cmod.reactions:
            if r.getUpperBound() != r.getLowerBound():
                r.setUpperBound(min(1000, r.getUpperBound()))
                r.setLowerBound(0)
        
        r = self.cmod.getReaction('R_Biomass_Ecoli_core_N_LPAREN_w_FSLASH_GAM_RPAREN__Nmet2')
        print(r.getUpperBound())
        print(r.getLowerBound())
        r.setLowerBound(0.8)
        r.setUpperBound(0.8)
                
        self.cmod.buildStoichMatrix()
        
        matrix = Matrix()
        labels, _ = matrix.addMetabolicNetwork(self.cmod)
        self.matroid = matroid.fromMatrix(matrix, labels)
        
        self.decomp = Decomposition(self.matroid)
        self.decomp.splitSimple()
        for v in self.decomp.listNonLeaves():
            if not v.isSimple():
                sim = v.matroid.similarityMatrix("combined")
                self.decomp.poolmanMethod(v, sim)

        write(self.decomp.makeGraphViz(cmod=self.cmod),"ecoli_textbook_bigg_split_poolman_combined.gv")
        #self.assertTrue(self.decomp.verifyFullyBranched())
        #self.assertTrue(self.decomp.verifyEdgeWidth())
        print("width = %d" % self.decomp.getWidth())

    def tearDown(self):
        pass

    def testEnumerateSubnetwork(self):
        CBSolver.FluxVariabilityAnalysis(self.cmod, objF2constr=False)
        v = self.decomp.getVertex(247)
        exclude = self.decomp.getVertex(283)
        enum = EFMEnumerator(self.cmod, v, exclude)
        F = enum.enumerateMinimal()
        print(F)
        lps = self.buildStrictLP(enum)
        #lpw = self.buildWeakLP(enum)
        modmatroid = self.matroid.contract([r for r in self.matroid.elems if r not in enum.module and r not in enum.fixed])
        #print(modmatroid.basis)
        #print(modmatroid.nbasis)
        print(enum.fixed)
        res = []
        for f in F:
            if self.isVertex(modmatroid, f):
                if not enum.isFeasible(f):
                    assert False
                sol = enum.lp.getSolution()
                interface = None
                for r in enum.module:
                    col = self.cmod.N.getColsByName(r)
                    if interface == None:
                        interface = col * sol[r]
                    else:
                        interface += col * sol[r]
                        
                # check strict feasibility
                res.append( (interface, self.isFeasible(enum, lps, f)) )
                print(f)
                for m in range(len(interface)):
                    v = interface[m]
                    if abs(v) > 1e-4:
                        print('%s: %g' % (self.cmod.N.row[m], v))
        #print(res)
            
    def isVertex(self, modmatroid, face):
        return modmatroid.isIndependent(face)
    
    def isFeasible(self, enum, lp, face):
        bounds = {}
        for r in enum.module:
            if r in face:
                # set lower bound to self.minflux
                lb = max(1e-4, self.cmod.getReactionLowerBound(r))
                ub = self.cmod.getReactionUpperBound(r)
                if lb > ub:
                    return False
                bounds[r] = (lb, ub)
            else:
                # fix flux to zero
                lb = self.cmod.getReactionLowerBound(r)
                if lb > 0:
                    return False
                bounds[r] = (lb, 0)
 
        lp.setBounds(bounds)
        lp.solve()
        
        if lp.getSolutionStatus() != 'LPS_UNDEF':
            return lp.isFeasible()
        else:
            print('warning: unable to solve LP')
            return True # if we keep it, we don't miss solutions
    
    def buildWeakLP(self, enum):
        lp = CBSolver.createSolver(self.cmod)
        for r in self.mnet.reactions:
            if r.getId() not in self.module:
                # remove all bounds from the variable
                # unless it has a fixed flux rate
                span = r.getFVAdata()[3]
                if (span > self.tol or isnan(span)):
                    # remove flux bounds
                    lp.setLowerBounds({r.getId():None})
                    lp.setUpperBounds({r.getId():None})
        # clear objective function
        lp.setObjective(reset=True)
        return lp

    def buildStrictLP(self, enum):
        lp = CBSolver.createSolver(self.cmod)
        # clear objective function
        lp.setObjective(reset=True)
        return lp

    @unittest.skip("takes quite long")
    def testEnumerate(self):
        CBSolver.FluxVariabilityAnalysis(self.cmod, objF2constr=False)
        enum = EFMEnumerator(self.cmod, self.decomp.root, None)
        fObj = CBModel.FluxObjective('R_EX_glc_e_rev_obj', reaction='R_EX_glc_LPAREN_e_RPAREN__rev', coefficient=1)
        objF = CBModel.Objective('minimize_nutrient', 'min')
        objF.addFluxObjective(fObj)
        self.cmod.addObjective(objF, active=True)
        CBSolver.analyzeModel(self.cmod)
        #CBSolver.MinimizeSumOfAbsFluxes(self.cmod, objF2constr=False)
        face = set()
        for r in self.cmod.reactions:
            if r.getValue() > 1e-5:
                face.add(r.getId())
                print('%s: %g' % (r.getId(), r.getValue()) )
        self.assertTrue(enum.test(face))
        F = enum.enumerateMinimal()
        print(F)

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