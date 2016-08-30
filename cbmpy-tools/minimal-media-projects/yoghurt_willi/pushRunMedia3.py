import os, numpy, time, itertools
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import cbmpy as cbm
import pprint

PPR = pprint.PrettyPrinter()

# I assume that the models live in a subdirectory of the directory where this script is run from called 'models'
modelDir = os.path.join(cDir, 'models')
constraintDir = os.path.join(cDir, 'constraints')
minimizationDir = os.path.join(cDir, 'minimizationsets')
# we assume there is only one model in 'models' in any case we use the first
#modelFile = [f for f in os.listdir(modelDir) if f.endswith('.xml')][0]

#tempoutF = 'd:\\@Dropbox\\Dropbox\\work\\temp_resq'

# bgoli-130813 for present purposes we hardcode 2 input files
#modelFile = '130726BP_NoTCA_split.xml'
modelFile = 'bulgTherSrcOnlyBulg.xml'

# just what it says whether to load an L3 file
USE_SBML3 = True

# bgoli-130813 this sets the roundoff factor for the MILP which has a tolereance of 1.0e-6 (default)
milpRoundOff = 8

# FVA parameters
FVA_TOL = None
FVA_RoundOff = 8

# this enables using a ranking strategy when testing possible media combinations smallest --> largest span
ENABLE_FVA_RANK = True

# magic perturbation constant
MINVAL = 1.0e-2
REVERSE_MINVAL = True

# get rid of pseudo alternatives
ENABLE_MAGIC_RESULT_FILTER = False

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
# 20130628 fixed Mac unicode insanity
sbmlF = str(os.path.join(cDir, 'models', modelFile))

if USE_SBML3:
    cmodBase = cbm.CBRead.readSBML3FBC(sbmlF)
else:
    cmodBase = cbm.CBRead.readSBML2FBA(sbmlF)


cbm.CBTools.processSBMLAnnotationNotes(cmodBase, annotation_key='note')
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


import copy
# GLOBALS
RESQ = []
RIDCACHE = []
CNTR = 0
CACHE = {}
DEBF = file('Debug-({}).txt'.format(modelFile), 'w')

def RUN(cmod3, minC, ignore, objFuncR, cntr):
    global RESQ, MINVAL, REVERSE_MINVAL, CACHE, DEBF, ENABLE_FVA_RANK
    cntr += 1
    ignore = copy.deepcopy(ignore)
    cbm.CBSolver.analyzeModel(cmod3, build_n=True)
    minCnt0 = 999999
    if cmod3.SOLUTION_STATUS != 'LPS_OPT':
        print('fba exit')
        DEBF.write('{}, infeasible model ({})\n'.format(cntr, cmod3.SOLUTION_STATUS))
        return {}
    try:
        minCnt0, lpx = cbm.CBSolver.cplx_MinimizeNumActiveFluxes(cmod3, selected_reactions=objFuncR, return_lp_obj=True,\
                                                                   quiet=True, objF2constr=False, pre_opt=True)
        if lpx.solution.get_status() != lpx.solution.status.MIP_optimal:
            print('mincnt fail')
            DEBF.write('{}, mincount not feasible ({})\n'.format(cntr, lpx.solution.get_status()))
            del lpx
            raise RuntimeError
        del lpx
    except:
        print('mincnt fail')
        DEBF.write('{}, mincount solve fail\n'.format(cntr))
        return {}

    if minCnt0 > minC:
        print('mincnt exit')
        DEBF.write('{}, mincount not min\n'.format(cntr))
        return {}

    nonZero = []
    for r_ in objFuncR:
        if abs(round(cmod3.getReaction(r_).getValue(), milpRoundOff)) > 0.0:
            nonZero.append(r_)

    print('nonZeroCnt: {}'.format(len(nonZero)))

    ## if len(nonZero) == 0:
        ## print('nonZerocnt exit')
        ## DEBF.write('{}, nonZeroCount is zero\n'.format(cntr))
        ## return {}
        

    testable = []
    for r_ in nonZero:
        if r_ not in ignore:
            testable.append(r_)
    #print('filipeCnt: {}'.format(len(testable)))
    print('Testable reactions:\n{}\n'.format(testable))

    if ENABLE_FVA_RANK:
        fva, fvan = cbm.CBSolver.FluxVariabilityAnalysis(cmod3, selected_reactions=testable, pre_opt=False,\
                                                tol=FVA_TOL, objF2constr=False, rhs_sense='lower', optPercentage=100.0,\
                                                quiet=True, markupmodel=False)
        testable_ranked = []
        for rr_ in testable:
            testable_ranked.append((fva[fvan.index(rr_)][4], fvan.index(rr_)))

        testable_ranked.sort()
        testable_ranked = [testable[a[1]] for a in testable_ranked]
        print('Testable reactions (ranked):\n{}\n'.format(testable_ranked))
        if testable != testable_ranked:
            DEBF.write('{},\"{}\",\"{}\"\n'.format(cntr, testable, testable_ranked))
        testable = testable_ranked

    res = {}
    if len(testable) == 0:
        DEBF.write('{}, no testable\n'.format(cntr))
        return {}
    else:
        for r_ in testable:
            cmod4 = cmod3.clone()
            cmod4.setReactionBounds(r_, 0.0, 0.0)
            for rr_ in testable:
                if rr_ != r_:
                    #rval = cmod4.getReaction(rr_).getValue()
                    if rr_ in neg_input_names or REVERSE_MINVAL:
                        DEBF.write('{}, setting upper bound ({}): {}\n'.format(cntr, rr_, -MINVAL))
                        cmod4.setReactionUpperBound(rr_, -MINVAL)
                    else:
                        DEBF.write('{}, setting lower bound ({}): {}\n'.format(cntr, rr_, MINVAL))
                        cmod4.setReactionLowerBound(rr_, MINVAL)
                    ignore.append(rr_)
            res['{}'.format(r_)] = RUN(cmod4, minC, ignore, objFuncR, cntr)
            #res['{}({})'.format(r_, cntr)] = RUN(cmod4, minC, ignore, objFuncR, cntr)
    del cmod3, cmod4
    print('Returning sub-iteration: {}'.format(cntr))
    DEBF.write('{}, success\n'.format(cntr))
    return res

xdata = {}
time0 = time.time()
for cset in constraintFiles:
    xdata[cset] = {}
    # first we clone the base model
    cmod = cmodBase.clone()
    # read bounds from file
    newBnds = readBounds(os.path.join(constraintDir, cset), delimiter=',')
    constraint_bound_reactions = []
    # set bounds from constraint file

    for b_ in newBnds:
        cmod.setReactionBounds(b_[0], b_[1], b_[2])
        constraint_bound_reactions.append(b_[0])

    # solve new model and make sure that it has a solution
    cbm.CBSolver.analyzeModel(cmod)

    # bgoli-130813 changed to not fail if there is a solver fail simply write empty result file
    #assert cmod.SOLUTION_STATUS == 'LPS_OPT', '\nSolver fail for cset: %s' % cset
    if cmod.SOLUTION_STATUS == 'LPS_OPT':
        # now we clone the new modified model, paranoia and because we can
        #efficiency# cmod = cmod.clone()
        cmod2 = cmod.clone()
        # new minset option
        for mset in minimizationSets:
            xdata[cset][mset] = {}
            objFuncReactions = []
            objFuncReactions = readMinimizationSet(os.path.join(minimizationDir, mset))


            # test for palsson back-to-front reactions in the minimization set
            NEG_INPUT = False
            neg_input_names = []
            for r_ in objFuncReactions:
                b_ = cmod.getReactionBounds(r_)
                if b_[1] < 0.0 and b_[2] <= 0.0:
                    NEG_INPUT = True
                    print('Negative input detected: {}'.format(b_[0]))
                    neg_input_names.append(b_[0])
            if NEG_INPUT:
                print('\n\nNegative input detected, for reactions:')
                print(neg_input_names)
                print('\n\nswitching to reversible (negative mode) ... ')
                ## print('... exiting, please reformulate the problem\n')
                ## os.sys.exit(-1)
            del NEG_INPUT
            neg_input_names = tuple(neg_input_names)
                
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
            P_maxpop = 1000
            P_timeout = 3600

            minCNt00 = cbm.CBSolver.cplx_MinimizeNumActiveFluxes(cmod, selected_reactions=objFuncReactions,\
                                                                               quiet=True, objF2constr=False, pre_opt=False)
            ignore = []
            cntr = 0


            fva, fvan = cbm.CBSolver.FluxVariabilityAnalysis(cmod, selected_reactions=objFuncReactions, pre_opt=False,\
                                                    tol=FVA_TOL, objF2constr=False, rhs_sense='lower', optPercentage=100.0,\
                                                    quiet=True, markupmodel=False)
            cbm.CBWrite.writeFVAtoCSV(fva, fvan, 'fva_({})_({})'.format(modelFile, cset), fbaObj=cmod)
            #blocked = []
            required = []
            testable = []
            for r_ in range(len(fvan)):
                if fvan[r_] in objFuncReactions:
                    #if round(fva[r_, 2], milpRoundOff) == round(fva[r_,3], milpRoundOff) and round(fva[r_, 2], milpRoundOff) <= 1.0e-10:
                    ##if abs(round(fva[r_, 2], milpRoundOff)) == abs(round(fva[r_,3], milpRoundOff)):
                        #blocked.append(fvan[r_])
                    # this switches to negative mode
                    if fvan[r_] in neg_input_names:
                        print(fva[r_])
                        if round(fva[r_, 3], FVA_RoundOff) < 0.0:
                            required.append(fvan[r_])
                        else:
                            testable.append(fvan[r_])
                    else:
                        if round(fva[r_, 2], FVA_RoundOff) > 0.0:
                            required.append(fvan[r_])
                        else:
                            testable.append(fvan[r_])

            DEBF.write('({})-({})\n'.format(modelFile, cset))

            print(neg_input_names)
            #print(objFuncReactions)
            print(required)
            print(testable)
            print(len(objFuncReactions), len(required), len(testable))
            time.sleep(5)
            
            #run_res = RUN(cmod2, minCNt00, required, objFuncReactions, cntr)
            run_res = RUN(cmod2, minCNt00, required, testable, cntr)

            PPR.pprint(run_res)
            magic_filtered = []
            if ENABLE_MAGIC_RESULT_FILTER:
                print('\nINFO: Magic result filter enabled!')
                for rk_ in list(run_res.keys()):
                    if len(run_res[rk_]) < 1:
                        run_res.pop(rk_)
                        magic_filtered.append(rk_)
                if len(run_res) < 1:
                    DEBF.write('No results returned from search algorithm, media results unreliable and generated for debugging purpose only!\n')
                    DEBF.write('For example, check that biomass constraints are not set too tight :->\n')
            PPR.pprint(run_res)

            xdata[cset][mset]['media'] = run_res
            xdata[cset][mset]['required'] = required
            xdata[cset][mset]['bounds'] = newBnds
            xdata[cset][mset]['model'] = cmod
            xdata[cset][mset]['objreact'] = objFuncReactions
            xdata[cset][mset]['magic_filtered'] = magic_filtered
            print('\nFiltered results: {}\n'.format(magic_filtered))
            time.sleep(3)
    else:
        print('\n\nSerious error 1 in constraint set: {}\n'.format(cset))
        os.sys.exit(-1)

DEBF.flush()
DEBF.close()

PPR.pprint(xdata)
F = file('media_search_results-({}).txt'.format(modelFile), 'w')
F.write(PPR.pformat(xdata))
F.close()

print('\nTotal results: {}\nTime: {} mins\n\n'.format(len(RESQ), (time.time()-time0)/60.))

#os.sys.exit(-2)

def unwrap(d):
    s = ''
    for k in d:
        if len(d[k]) == 0:
            s += '{},'.format(k)
        else:
            s += '{},{}'.format(k, unwrap(d[k]))
    return s

lines = []

for c in xdata:
    L = 'MODEL,{},'.format(modelFile)
    lines.append(L[:-1])
    L = 'CONSTRAINTS,{},'.format(c)
    lines.append(L[:-1])

    for m in xdata[c]:
        L = 'MINSET,{},'.format(m)
        lines.append(L[:-1])
        L2 = 'REQUIRED,'
        for r in xdata[c][m]['required']:
            L2 += '{},'.format(r)
        lines.append(L2[:-1])
        if 'magic_filtered' in xdata[c][m]:
            L2a = 'FILTERED,'
            for r in xdata[c][m]['magic_filtered']:
                L2a += '{},'.format(r)
        lines.append(L2a[:-1])
        cntr = 1
        groups = []
        for r in xdata[c][m]['media']:
            L3 = 'ALT{},{},'.format(cntr,r)
            cntr += 1
            S = unwrap(xdata[c][m]['media'][r])
            L3 += '{}'.format(S)
            lines.append(L3[:-1])
            g2 = [s.strip() for s in L3[:-1].split(',')]
            g2.pop(0)
            groups.append(g2)
        combis = list(itertools.product(*groups))
        xdata[c][m]['combis'] = combis
        print('Media combinations ({})({}): {}'.format(c, m, len(combis)))

cbm.CBTools.storeObj(xdata, 'media_search_results-({})'.format(modelFile))

try:
    F = file('media_search_results-({}).csv'.format(modelFile), 'w')
    for l in lines:
        F.write('{}\n'.format(l))
    F.close()
except:
    print('Close Excel file')


mediaList = xdata.keys()
objList = xdata[mediaList[0]].keys()

print(mediaList)
print(objList)

# ['media', 'model', 'required', 'bounds']
ydata = {}
sdata = {}
ydata[modelFile] = {}
sdata[modelFile] = {}
for o_ in objList:
    ydata[modelFile][o_] = {}
    sdata[modelFile][o_] = {}
    for m_ in mediaList:
        ydata[modelFile][o_][m_] = {}
        sdata[modelFile][o_][m_] = {}
        try:
            cmod = xdata[m_][o_]['model']
            #cmod.setReactionLowerBound('R_BP_biomass_130704', 0.0)
            #cmod.setReactionUpperBound('R_BP_biomass_130704', float('inf'))
            newBnds = xdata[m_][o_]['bounds']
            objreact = xdata[m_][o_]['objreact']
            required = xdata[m_][o_]['required']
            combis = xdata[m_][o_]['combis']
            for b_ in objreact:
                cmod.setReactionBounds(b_, 0.0, 0.0)
            for b_ in required:
                if b_ in neg_input_names or REVERSE_MINVAL:
                    cmod.setReactionBounds(b_, -float('inf'), 0.0)
                else:
                    cmod.setReactionBounds(b_, 0.0, float('inf'))
            mcntr = 1
            for c_ in combis:
                for cb_ in c_:
                    if cb_ in neg_input_names or REVERSE_MINVAL:
                        cmod.setReactionBounds(cb_, -float('inf'), 0.0)
                    else:
                        cmod.setReactionBounds(cb_, 0.0, float('inf'))
                cbm.CBCPLEX.cplx_MinimizeSumOfAbsFluxes(cmod, objF2constr=False, pre_opt=False, quiet=True)
                media = {}
                solution = {}
                for cb_ in c_:
                    media[cb_] = cmod.getReaction(cb_).getValue()
                    cmod.setReactionBounds(cb_, 0.0, 0.0)
                for rb_ in required:
                    media[rb_] = cmod.getReaction(rb_).getValue()
                rids = cmod.getReactionIds()
                rids.sort()
                for r_ in rids:
                    solution[r_] = cmod.getReaction(r_).getValue()
                ydata[modelFile][o_][m_][mcntr] = media
                sdata[modelFile][o_][m_][mcntr] = solution
                mcntr += 1
        except Exception as ex:
            print(ex)
            print('\nMedia generation failed for: {} : {}'.format(o_, m_))
            #os.sys.exit()

PPR.pprint(ydata)
F = file('media_final-({}).txt'.format(modelFile), 'w')
F.write(PPR.pformat(ydata))
F.close()
F = file('solutions-({}).txt'.format(modelFile), 'w')
F.write(PPR.pformat(sdata))
F.close()

cbm.CBTools.storeObj(ydata, 'media_final-({})'.format(modelFile))
cbm.CBTools.storeObj(sdata, 'solutions-({})'.format(modelFile))

for m_ in ydata:
    for o_ in ydata[m_]:
        for c_ in ydata[m_][o_]:
            mset = ydata[m_][o_][c_]
            keys = mset.keys()
            print keys
            keys.sort()
            if len(keys) >= 1:
                F = file('Media_({})_({}_({})_({})).csv'.format(m_, o_, c_, max(keys)), 'w')
                F.write('{}\n'.format(c_))
                for k_ in keys:
                    for r_ in mset[k_]:
                        F.write('{},{},{}\n'.format(k_, r_, mset[k_][r_]))
            else:
                F = file('Media_({})_({}_({})_({})).csv'.format(m_, o_, c_, 0), 'w')
            F.close()

for m_ in sdata:
    for o_ in sdata[m_]:
        for c_ in sdata[m_][o_]:
            mset = sdata[m_][o_][c_]
            keys = mset.keys()
            print keys
            keys.sort()
            if len(keys) >= 1:
                out = []
                for k_ in keys:
                    out.append([])
                    for r_ in mset[k_]:
                        out[-1].append((r_, mset[k_][r_]))
                sdata[m_][o_][c_] = out
                F = file('Solutions_({})_({}_({})_({})).csv'.format(m_, o_, c_, max(keys)), 'w')
                L = '{},'.format(c_)
                for f_ in range(len(out)):
                    L += '{},'.format(f_+1)
                F.write('{}\n'.format(L[:-1]))
                try:
                    for r_ in range(len(out[0])):
                        L = '{},'.format(out[0][r_][0])
                        for d_ in range(len(out)):
                            L += '{},'.format(out[d_][r_][1])
                        F.write('{}\n'.format(L[:-1]))
                except:
                    pass
            else:
                F = file('Solutions_({})_({}_({})_({})).csv'.format(m_, o_, c_, 0), 'w')
            F.close()

cbm.CBTools.storeObj(sdata, 'solutions-({})'.format(modelFile))

endtime = time.time()
F = file('TimeCustomSearch-({}).txt'.format(modelFile), 'w')
F.write('Total analysis time (seconds): {}\n'.format(round(endtime-time0,1)))
F.close()
print('\n\nTotal analysis time (seconds): {}'.format(round(endtime-time0,1)))
