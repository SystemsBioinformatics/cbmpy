import os, numpy, time
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import cbmpy as cbm

# I assume that the models live in a subdirectory of the directory where this script is run from called 'models'
modelDir = os.path.join(cDir, 'models')
constraintDir = os.path.join(cDir, 'constraints')
minimizationDir = os.path.join(cDir, 'minimizationsets')
# we assume there is only one model in 'models' in any case we use the first
modelFile = [f for f in os.listdir(modelDir) if f.endswith('.xml')][0]

# bgoli-130813 for present purposes we hardcode 2 input files
#modelFile = '130726BP_NoTCA_split.xml'
#modelFile = '130726BP_TCA_split.xml'
modelFile = 'bulgTherSrcOnlyBulg.xml'

# set the SBML level of the input file  (True=3)
USE_SBML3 = True

# bgoli-130813 this sets the roundoff factor for the MILP which has a tolereance of 1.0e-6 (default)
milpRoundOff = 8
# bgoli-130813 use "suboptimal" biomass [biomassHackPC*LB : biomassHackPC*UB] None means disabled
biomassHack = None

# multiple constraints possible
constraintFiles = []
for f in os.listdir(constraintDir):
    if f.endswith('.csv'):
        constraintFiles.append(f.strip())
assert len(constraintFiles) > 0, '\nNo constraint sets defined'

# multiple minimization sets
minimizationSets = []
for f in os.listdir(minimizationDir):
    if f.endswith('.csv'):
        minimizationSets.append(f.strip())
assert len(minimizationSets) > 0, '\nNo minimization sets defined'

#constraintFiles = ['C_source_screen.csv']

print '\nmodelFile'
print modelFile
print '\nconstraintFiles'
print constraintFiles

# ok now lets set up our model

# let's set up the base model
if USE_SBML3:
    cmodBase = cbm.CBRead.readSBML3FBC(modelFile, os.path.join(cDir, 'models'))
else:
    cmodBase = cbm.CBRead.readSBML2FBA(modelFile, os.path.join(cDir, 'models'))
#cbm.CBTools.processSBMLAnnotationNotes(cmodBase, annotation_key='note')
cbm.CBSolver.analyzeModel(cmodBase)

# here we define a function that will read the constraint files
def readBounds(filename, delimiter=','):
    """
    Give the name of a bounds file and it returns a list of tuples: (Rid,lb,ub)

    """
    F = file(filename, 'r')
    bounds = []
    for L in F:
        l = L.split(delimiter)
        #print l
        if l[0].strip().replace('"','') == 'Reaction ID':
            pass
        elif len(l) == 3:
            bounds.append((l[0].strip().replace('"',''), float(l[1].strip()), float(l[2].strip())))
        else:
            raise RuntimeError, '\nInvalid bounds file format'
    F.close()
    assert len(bounds) > 0, '\nBound file read error, check line endings they should be windows!\n'
    return bounds

def readMinimizationSet(filename):
    """
    Give the name of a minimization file it returns a list of rids's: ['R1','R2',...]

    """
    F = file(filename, 'r')
    minset = []
    for L in F:
        Lstrip = L.strip()
        if 'Reaction ID' in Lstrip or Lstrip == '' or Lstrip == ' ':
            pass
        else:
            minset.append(Lstrip)

    F.close()
    assert len(minset) > 0, '\nBound file read error, check line endings\n'
    if len(minset) == 1:
        print('\nIs this the objective function set you expect: {}?'.format(minset))
    return minset

sdata = []

time0 = time.time()
for cset in constraintFiles:
    # first we clone the base model
    cmod = cmodBase.clone()
    # read bounds from file
    newBnds = readBounds(os.path.join(constraintDir, cset), delimiter=',')
    constraint_bound_reactions = []
    # set bounds from constraint file

    for b_ in newBnds:
        assert cmod.getReaction(b_[0]) != None, '\nERROR: %s is not a valid reaction id (%s)!' % (b_[0], cset)
        if biomassHack != None and b_[0] == 'R_BP_biomass_130704':
            print('Biomasshack enabled')
            cmod.setReactionBounds(b_[0], biomassHack*b_[1], biomassHack*b_[2])
            print(cmod.getReactionBounds(b_[0]))
            time.sleep(2)
        else:
            cmod.setReactionBounds(b_[0], b_[1], b_[2])
        constraint_bound_reactions.append(b_[0])

    # solve new model and make sure that it has a solution
    cbm.CBSolver.analyzeModel(cmod)

    # bgoli-130813 changed to not fail if there is a solver fail simply write empty result file
    #assert cmod.SOLUTION_STATUS == 'LPS_OPT', '\nSolver fail for cset: %s' % cset
    if cmod.SOLUTION_STATUS == 'LPS_OPT':
        # now we clone the new modified model, paranoia and because we canS
        cmod2 = cmod.clone()

        # new minset option
        for mset in minimizationSets:
            objFuncReactions = []
            objFuncReactions = readMinimizationSet(os.path.join(minimizationDir, mset))

            minSum = minCnt = 0.0

            # now we do a flux sum minimization on cmod
            #minSum = cbm.CBSolver.cplx_MinimizeSumOfAbsFluxes(cmod, selected_reactions=objFuncReactions, quiet=True)
            minSum = cbm.CBSolver.cplx_MinimizeSumOfAbsFluxes(cmod, selected_reactions=objFuncReactions, objF2constr=False, pre_opt=False, quiet=True)


            # now lets do a minimum number active flux minimization

            # bgoli-130813 this sets the relative gap to a solution, normally should be 0.0
            P_relgap = 0.0
            #P_relgap = 1.0e-3

            # bgoli-130813 this sets a maximum of 1 hour per solve and a maximum of 1000 solutions
            #P_maxpop = 100
            #P_timeout = 300
            P_maxpop = 2100000000
            P_timeout = 3600

            popSol, minCnt, lpx = cbm.CBSolver.cplx_MinimizeNumActiveFluxes(cmod2, selected_reactions=objFuncReactions,\
                                                               quiet=True, populate=(P_relgap,P_maxpop, P_timeout),\
                                                               pre_opt=False, objF2constr=False, return_lp_obj=True)
            ## Debug for brett
            lpDir = os.path.join(cDir, 'lpx')
            if not os.path.exists(lpDir):
                os.makedirs(lpDir)
            lpf = os.path.join(lpDir, '(%s)-(%s)-(%s)-(%s).lp' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt))
            lpx.write(lpf)
            lpf = os.path.join(lpDir, '(%s)-(%s)-(%s)-(%s).xml' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt))
            if USE_SBML3:
                cbm.CBWrite.writeSBML3FBC(cmod, lpf)
            else:
                cbm.CBWrite.writeSBML2FBA(cmod, lpf)
            ## Debug for brett

            print('\n***\nMinimum sum of fluxes is: {}\nMinimum active flux count is: {}\nNumber of solutions: {}\n'.format(minSum, minCnt, len(popSol)-1))
            time.sleep(2)

            # we can now start assembling data
            minSumData = []
            baseData = []
            minNumData = []
            for r_ in cmod.getReactionIds():
                baseData.append(cmod.getReaction(r_).getValue())
                minSumData.append(cmod.getReaction(r_).getValue())
                ridx = popSol[0].index(r_)
                rpopData = []
                for s_ in range(1, len(popSol)):
                    rpopData.append(popSol[s_][ridx])
                minNumData.append(rpopData)

            #print objFuncReactions
            #print baseData
            #print minSumData
            #print minNumData
            #print cmod.getActiveObjective().getValue()
            #print cmod2.getActiveObjective().getValue()


            # and write it all out, what we are aiming for is
            # J,OptFlux,MinSumFlux,MinActive1,MinActive2,MinActive3 ...

            if biomassHack != None:
                if P_relgap != 0.0:
                    F = file(os.path.join(cDir, 'Result_meth0-(%s)-(%s)-(%s)-(%s)-(%s)-(biomasshack=%s).csv' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt, len(popSol)-1, biomassHack, P_relgap)),'w')
                else:
                    F = file(os.path.join(cDir, 'Result_meth0-(%s)-(%s)-(%s)-(%s)-(%s)-(biomasshack=%s)-(P_relgap=%s).csv' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt, len(popSol)-1, biomassHack)),'w')
            else:
                if P_relgap != 0.0:
                    F = file(os.path.join(cDir, 'Result_meth0-(%s)-(%s)-(%s)-(%s)-(%s)-(P_relgap=%s).csv' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt, len(popSol)-1, P_relgap)),'w')
                else:
                    F = file(os.path.join(cDir, 'Result_meth0-(%s)-(%s)-(%s)-(%s)-(%s).csv' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt, len(popSol)-1)),'w')
            head = 'R,OptFlux,MinSumFlux,Constraint,Minset'
            for s_ in range(1, len(popSol)):
                head += ',MinActFlux%i' % s_
            F.write('%s\n' % head)

            reactionList2 = cmod.getReactionIds()
            for r_ in range(len(reactionList2)):
                line = reactionList2[r_]
                line += ',%s' % baseData[r_]
                line += ',%s' % minSumData[r_]
                if reactionList2[r_] in constraint_bound_reactions:
                    line += ',1'
                else:
                    line += ',0'
                if reactionList2[r_] in objFuncReactions:
                    line += ',1'
                else:
                    line += ',0'
                for d_ in minNumData[r_]:
                    line += ',%s' % round(d_, milpRoundOff)
                F.write('%s\n' % line)
            F.close()
            sdata.append((baseData, minSumData, minNumData, minSum, minCnt, popSol, (cset, mset)))
    else:
        for mset in minimizationSets:
            minCnt = 0
            F = file(os.path.join(cDir, 'Result_meth0-(%s)-(%s)-(%s)-(%s)-(SOLVERFAIL).csv' % (modelFile.replace('.xml',''), cset.replace('.csv',''), mset.replace('.csv',''), minCnt)),'w')
            F.close()
    del cmod, cmod2

cbm.CBTools.storeObj(sdata, 'rm_data-({})'.format(modelFile))

endtime = time.time()
F = file('TimeOptimizedMILP-({}).txt'.format(modelFile), 'w')
F.write('Total analysis time (seconds): {}\n'.format(round(endtime-time0,1)))
F.close()
print('\n\nTotal analysis time (seconds): {}'.format(round(endtime-time0,1)))

















"""
#for model in models:
    #time0 = time.time()
    #cmod = cbm.CBRead.readSBML2FBA(model, modelDir)
    #cbm.CBTools.processSBMLAnnotationNotes(cmod)
    #cbm.CBSolver.analyzeModel(cmod, oldlpgen=True)
    #if not useMultiCore:
        #fvao, fvaon = cbm.CBSolver.FluxVariabilityAnalysis(cmod, selected_reactions=objFuncReactions, optPercentage=100.0, quiet=True, oldlpgen=True)
    #else:
        #fvao, fvaon = cbm.CBMultiCore.runMultiCoreFVA(cmod, selected_reactions=objFuncReactions, optPercentage=100.0, procs=coreNumber)
    ## results are written to a file <modelname>-(<optimum>).csv
    #cbm.CBWrite.WriteFVAtoCSV(model.replace('.xml', '-(%s)' % round(cmod.getActiveObjective().getValue(), 3)), fvao, fvaon, fbaObj=cmod)
    #print '\nModel %s FVA took %s seconds' % (model, round(time.time()-time0, 2))
"""


