"""
CBMPy: CBGLPK module
====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2015 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBGLPK.py 390 2015-10-05 13:44:45Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

HAVE_SYMPY = False
try:
    import sympy
    if int(sympy.__version__.split('.')[1]) >= 7 and int(sympy.__version__.split('.')[2]) >= 4:
        HAVE_SYMPY = True
    else:
        del sympy
except ImportError:
    HAVE_SYMPY = False

import os, time
import numpy
from . import CBWrite
from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


HAVE_GLPK = False
GLPK_SOLUTION_STATUS = None

try:
    import glpk
    lp = glpk.LPX()
    HAVE_GLPK = True
except:
    raise ImportError

# configuration options for GLPK
GLPK_CFG = {'simplex' : {'meth' : glpk.LPX.PRIMAL,
                       #'meth' : glpk.LPX.DUAL,
                       #'meth' : glpk.LPX.DUALP,
                       'tol_bnd' : 1.0e-6,
                       'tol_dj'  : 1.0e-6,
                       'tol_piv' : 1.0e-10
                       }
          }

GLPK_STATUS = {
    1 : 'LPS_UNDEF',
    2 : 'LPS_FEAS',
    3 : 'LPS_INFEAS',
    4 : 'LPS_NOFEAS',
    5 : 'LPS_OPT',
    6 : 'LPS_UNBND'}

GLPK_STATUS2 = {
    'opt' : 'LPS_OPT',
    'undef' : 'LPS_UNDEF',
    'feas' : 'LPS_FEAS',
    'infeas' : 'LPS_INFEAS',
    'nofeas' : 'LPS_NOFEAS',
    'unbnd' : 'LPS_UNBND'}


GLPK_SILENT_MODE = True
GLPK_INFINITY = 1.0e9

def glpk_constructLPfromFBA(fba, fname=None):
    """
    Create a GLPK LP in memory.
    - *fba* an FBA object
    - *fname* optional filename if defined writes out the constructed lp

    """
    _Stime = time.time()

    # define model and add variables
    lp = glpk.LPX()
    lp.name = fba.getPid()
    lp.cols.add(len(fba.N.col))


    if HAVE_SYMPY and fba.N.__array_type__ == sympy.MutableDenseMatrix:
        print('INFO: GLPK requires floating point, converting N')
        Nmat = numpy.array(fba.N.array).astype('float')
        RHSmat = numpy.array(fba.N.RHS).astype('float')
        if fba.CM != None:
            CMmat = numpy.array(fba.CM.array).astype('float')
            CMrhs = numpy.array(fba.CM.RHS).astype('float')
    else:
        Nmat = fba.N.array
        RHSmat = fba.N.RHS
        if fba.CM != None:
            CMmat = fba.CM.array
            CMrhs = fba.CM.RHS

    varMap = {}
    for n_ in range(Nmat.shape[1]):
        varMap[fba.N.col[n_]] = n_
        lp.cols[n_].name = fba.N.col[n_]

    #print varMap

    try:
        # define objective
        osense = fba.getActiveObjective().operation.lower()
        if osense in ['minimize', 'minimise', 'min']:
            lp.obj.maximize = False
        elif osense in ['maximize', 'maximise', 'max']:
            lp.obj.maximize = True
        else:
            raise RuntimeError('\n%s - is not a valid objective operation' % osense)
        lp.obj.name = fba.getActiveObjective().getPid()
        for fo_ in fba.getActiveObjective().fluxObjectives:
            lp.obj[varMap[fo_.reaction]] = fo_.coefficient
    except AttributeError:
        print('\nWARNING(GLPK create LP): no objective function defined')
    # create N constraints
    lp.rows.add(Nmat.shape[0])
    #conMap = {}
    #for n_ in range(Nmat.shape[0]):
        ##conMap[fba.N.row[n_]] = n_
        #lp.rows[n_].name = fba.N.row[n_]

    #tnew = time.time()
    for r_ in range(Nmat.shape[0]):
        # name and coefficients
        newCon = []
        for c_ in range(Nmat.shape[1]):
            newCon.append((c_, Nmat[r_,c_]))
        lp.rows[r_].name = fba.N.row[r_]
        lp.rows[r_].matrix = newCon

        # sense and rhs
        rhs = RHSmat[r_]
        if fba.N.operators[r_] in ['<=','<','L']:
            lp.rows[r_].bounds = None, rhs
        elif fba.N.operators[r_] in ['>=','>','G']:
            lp.rows[r_].bounds = rhs, None
        elif fba.N.operators[r_] in ['=','E']:
            lp.rows[r_].bounds = rhs
        else:
            raise RuntimeError('INFO: invalid operator: %s' % fba.N.operators[n])

    # add user defined constraints
    if fba.CM != None:
        baseRows = len(lp.rows)
        lp.rows.add(CMmat.shape[0])
        for r_ in range(CMmat.shape[0]):
            # name and coefficients
            newCon = []
            for c_ in range(CMmat.shape[1]):
                newCon.append((c_, CMmat[r_, c_]))
            lp.rows[baseRows+r_].name = fba.CM.row[r_]
            lp.rows[baseRows+r_].matrix = newCon

            # sense and rhs
            rhs = CMrhs[r_]
            if fba.CM.operators[r_] in ['<=','<','L']:
                lp.rows[baseRows+r_].bounds = None, rhs
            elif fba.CM.operators[r_] in ['>=','>','G']:
                lp.rows[baseRows+r_].bounds = rhs, None
            elif fba.CM.operators[r_] in ['=','E']:
                lp.rows[baseRows+r_].bounds = rhs
            else:
                raise RuntimeError('INFO: invalid operator: %s' % fba.N.operators[n])

    # add bounds
    for r_ in fba.reactions:
        lb = ub = None
        lb = fba.getReactionLowerBound(r_.getPid())
        ub = fba.getReactionUpperBound(r_.getPid())

        if lb in ['Infinity', 'inf', 'Inf', 'infinity']:
            lb = GLPK_INFINITY
        elif lb in ['-Infinity', '-inf', '-Inf', '-infinity', None]:
            lb = -GLPK_INFINITY
        elif numpy.isinf(lb):
            if lb < 0.0:
                lb = -GLPK_INFINITY
            else:
                lb = GLPK_INFINITY
        if ub in ['Infinity', 'inf', 'Inf', 'infinity', None]:
            ub = GLPK_INFINITY
        elif ub in ['-Infinity', '-inf', '-Inf', '-infinity']:
            ub = -GLPK_INFINITY
        elif numpy.isinf(ub):
            if ub < 0.0:
                ub = -GLPK_INFINITY
            else:
                ub = GLPK_INFINITY

        if ub != GLPK_INFINITY and lb != -GLPK_INFINITY and ub == lb:
            lp.cols[varMap[r_.getPid()]].bounds = lb
        elif ub != GLPK_INFINITY and lb != -GLPK_INFINITY:
            lp.cols[varMap[r_.getPid()]].bounds = lb, ub
        elif ub != GLPK_INFINITY:
            lp.cols[varMap[r_.getPid()]].bounds = None, ub
        elif lb != -GLPK_INFINITY:
            lp.cols[varMap[r_.getPid()]].bounds = lb, None
        else:
            lp.cols[varMap[r_.getPid()]].bounds = None

    print('\ngplk_constructLPfromFBA time: {}\n'.format(time.time() - _Stime))
    if fname != None:
        lp.write(cpxlp=fname+'.lp')
    return lp

def glpk_Solve(lp, method='s'):
    """
    Solve the LP and create a status attribute with the solution status

     - *method* [default='s'] 's' = simplex, 'i' = interior, 'e' = exact

    GLPK solver options can be set in the dictionary GLPK_CFG

    """

    if method == 'i':
        glpksol = lp.interior()
    elif method == 'e':
        glpksol = lp.exact()
    else:
        glpksol = lp.simplex(**GLPK_CFG['simplex'])


    global GLPK_SOLUTION_STATUS
    GLPK_SOLUTION_STATUS = glpk_getSolutionStatus(lp)

    #if status == 'LPS_UNDEF':
        #sd.presolve = False
        #status = glpk.glp_simplex(lp, sd)
    #if status == 'LPS_UNDEF':
        #print('\nINFO: Primal solver failure switching to dual-primal solver\n')
        #sd.presolve = False
        #sd.obj_ul = 1.0e10
        #sd.obj_ll = 1.0e10
        #sd.meth = glpk.GLP_DUALP
        #glpk.glp_simplex(lp, sd)
        #status = glpk_getSolutionStatus(lp)



    if GLPK_SOLUTION_STATUS in ['LPS_UNDEF', 'LPS_FEAS', 'LPS_INFEAS', 'LPS_NOFEAS', 'LPS_OPT', 'LPS_UNBND']:
        if not GLPK_SILENT_MODE:
            print('Solution status returned as: {}'.format(GLPK_SOLUTION_STATUS))
            print("Objective value = " , lp.obj.value)
        return GLPK_SOLUTION_STATUS
    else:
        print("INFO: No solution available ({})".format(GLPK_SOLUTION_STATUS))
        return None

def glpk_getSolutionStatus(lp):
    """
    Returns one of:

     - *LPS_OPT*: solution is optimal;
     - *LPS_FEAS*: solution is feasible;
     - *LPS_INFEAS*: solution is infeasible;
     - *LPS_NOFEAS*: problem has no feasible solution;
     - *LPS_UNBND*: problem has unbounded solution;
     - *LPS_UNDEF*: solution is undefined.

    """
    return GLPK_STATUS2[lp.status]

def glpk_analyzeModel(f, lpFname=None, return_lp_obj=False, with_reduced_costs='unscaled', with_sensitivity=False,\
                      del_intermediate=False, build_n=True, quiet=False, oldlpgen=False, method='s'):
    """
    Optimize a model and add the result of the optimization to the model object
    (e.g. `reaction.value`, `objectiveFunction.value`). The stoichiometric
    matrix is automatically generated. This is a common function available in all
    solver interfaces. By default returns the objective function value

     - *f* an instantiated PySCeSCBM model object
     - *lpFname* [default=None] the name of the intermediate LP file saved when this has a string value.
     - *return_lp_obj* [default=False] off by default when enabled it returns the PyGLPK LP object
     - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
       or anything else which is interpreted as 'None'. Scaled means s_rcost = (r.reduced_cost*rval)/obj_value
     - *with_sensitivity* [default=False] add solution sensitivity information (not yet implemented)
     - *del_intermediate* [default=False] delete the intermediary files after updating model object, useful for server applications
     - *build_n* [default=True] generate stoichiometry from the reaction network (reactions/reagents/species)
     - *quiet* [default=False] suppress glpk output
     - *method* [default='s'] select the GLPK solver method, see the GLPK documentation for details

       - 's': simplex
       - 'i': interior
       - 'e': exact

    """

    if build_n:
        f.buildStoichMatrix()
        #CBTools.addStoichToFBAModel(f)
    fid = f.id

    if with_reduced_costs == 'scaled':
        f.SCALED_REDUCED_COSTS = True
    elif with_reduced_costs == 'unscaled':
        f.SCALED_REDUCED_COSTS = False
    else:
        f.SCALED_REDUCED_COSTS = None

    if lpFname == None:
        flp = glpk_constructLPfromFBA(f, fname=None)
        f.id = '_glpktmp_.tmp'
    else:
        flp = glpk_constructLPfromFBA(f, fname=lpFname)
        fid = lpFname

    _Stime = time.time()

    print('\nglpk_analyzeModel FBA --> LP time: {}\n'.format(time.time() - _Stime))
    f.id = fid
    glpk_Solve(flp, method)

    glpk_setFBAsolutionToModel(f, flp, with_reduced_costs=with_reduced_costs)

    f.SOLUTION_STATUS =  glpk_getSolutionStatus(flp)
    # TODO: need to synchronise all solvers to same integer system
    if f.SOLUTION_STATUS == 'LPS_OPT':
        f.SOLUTION_STATUS_INT = 1
    else:
        f.SOLUTION_STATUS_INT = 999

    if oldlpgen and del_intermediate:
        os.remove(LPF)
    objv = f.getActiveObjective().getValue()
    print('\nanalyzeModel objective value: {}\n'.format(objv))
    if return_lp_obj:
        return flp
    else:
        del flp
        return objv

def glpk_setFBAsolutionToModel(fba, lp, with_reduced_costs='unscaled'):
    """
    Sets the FBA solution from a CPLEX solution to an FBA object

     - *fba* and fba object
     - *lp* a CPLEX LP object
     - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
       or anything else which is interpreted as None. Scaled is: s_rcost = (r.reduced_cost*rval)/obj_value

    """
    sol, objname, objval = glpk_getOptimalSolution(lp)
    if glpk_getSolutionStatus(lp) == 'LPS_OPT':
        fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = sol, objval
    else:
        fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = sol, numpy.NaN
    for r in fba.reactions:
        rid = r.getPid()
        if rid in sol:
            r.value = sol[rid]
        else:
            r.value = None
    scaled = False
    if with_reduced_costs == 'scaled':
        scaled = True
        with_reduced_costs = True
    elif with_reduced_costs == 'unscaled':
        scaled = False
        with_reduced_costs = True
    else:
        with_reduced_costs = False

    if objval != None and sol != {} and with_reduced_costs:
        RC = glpk_getReducedCosts(lp, scaled=scaled)
        setReducedCosts(fba, RC)
    else:
        setReducedCosts(fba, {})

def glpk_getOptimalSolution(c):
    """
    From a GLPK model extract a tuple of solution, ObjFuncName and ObjFuncVal
    """
    s_val = []
    s_name = []
    fba_sol = {}
    objf_name = None
    objf_val = None
    try:
        objf_name = c.obj.name
        objf_val = c.obj.value
        for n in c.cols:
            fba_sol.update({n.name : n.value})
    except Exception as ex:
        print(ex)
        print('WARNING: No solution to get')
        s_val = []
        s_name = []
        fba_sol = {}
        objf_val = None
    del s_name,s_val
    return fba_sol, objf_name, objf_val

def glpk_getReducedCosts(c, scaled=False):
    """
    Extract ReducedCosts from LP and return as a dictionary 'Rid' : reduced cost

     - *c* a GLPK LP object
     - *scaled* scale the reduced cost by the optimal flux value

    """
    s_name = []
    r_costs = []
    s_val = []

    for c_ in c.cols:
        s_name.append(c_.name)
        r_costs.append(c_.dual)
        s_val.append(c_.value)
    objf_val = c.obj.value
    output = {}
    for v in range(len(s_name)):
        if scaled:
            try:
                r_val = r_costs[v]*s_val[v]/objf_val
            except:
                r_val = float('nan')
        else:
            r_val = r_costs[v]
        output.update({s_name[v] : r_val})
    del s_name, r_costs, s_val
    return output

def glpk_FluxVariabilityAnalysis(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=True, debug=False, oldlpgen=False, markupmodel=True, default_on_fail=False, roundoff_span=10, method='s'):
    """
    Perform a flux variability analysis on an fba model:

     - *fba* an FBA model object
     - *selected reactions* [default=None] means use all reactions otherwise use the reactions listed here
     - *pre_opt* [default=True] attempt to presolve the FBA and report its results in the ouput, if this is disabled and *objF2constr* is True then the vid/value of the current active objective is used
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)
     - *work_dir* [default=None] the FVA working directory for temporary files default = cwd+fva
     - *debug* [default=False] if True write out all the intermediate FVA LP's into work_dir
     - *quiet* [default=False] if enabled supress CPLEX output
     - *objF2constr* [default=True] add the model objective function as a constraint using rhs_sense etc. If
       this is True with pre_opt=False then the id/value of the active objective is used to form the constraint
     - *markupmodel* [default=True] add the values returned by the fva to the reaction.fva_min and reaction.fva_max
     - *default_on_fail* [default=False] if *pre_opt* is enabled replace a failed minimum/maximum with the solution value
     - *roundoff_span* [default=10] number of digits is round off (not individual min/max values)
     - *method* [default='s'] select the GLPK solver method, see the GLPK documentation for details

       - 's': simplex
       - 'i': interior
       - 'e': exact


    Returns an array with columns Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus and a list containing the row names.

    """
    if work_dir == None:
        work_dir = os.getcwd()
    else:
        assert os.path.exists(work_dir), '\nWhat did you think would happen now!'
    if debug:
        debug_dir = os.path.join(work_dir,'DEBUG')
        if not os.path.exists(debug_dir):
            os.mkdir(debug_dir)
    s2time = time.time()
    # generate a presolution
    print('GLPK is using solver option: "{}"'.format(method))
    cpx = OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = glpk_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, method=method)

    # if required add the objective function as a constraint
    if objF2constr:
        glpk_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, tol, optPercentage)
    if debug:
        cpx.write(cpxlp=os.path.join(debug_dir, 'FVA_base.lp'))

    # do the FVA
    NUM_FLX = len(fba.reactions)
    print('Total number of reactions: {}'.format(NUM_FLX))
    if selected_reactions != None:
        rids = fba.getReactionIds()
        for r in selected_reactions:
            assert r in rids, "\n%s is not a valid reaction name" % r
    else:
        selected_reactions = fba.getReactionIds()
    NUM_FLX = len(selected_reactions)
    print('Number of user selected variables: {}'.format(NUM_FLX))
    try:
        OUTPUT_ARRAY = numpy.zeros((NUM_FLX, 7), numpy.double)
    except AttributeError:
        OUTPUT_ARRAY = numpy.zeros((NUM_FLX, 7))
    OUTPUT_NAMES = []
    cntr = 0
    tcnt = 0
    # this is a memory hack --> prevents solver going berserk
    mcntr = 0
    mcntr_cntrl = 3
    mps_filename = '_{}_.mps'.format(str(time.time()).split('.')[0])
    cpx.write(mps=mps_filename)
    for Ridx in range(NUM_FLX):
        R = selected_reactions[Ridx]
        OUTPUT_NAMES.append(R)
        max_error_iter = 1
        GOMIN = True
        gomin_cntr = 0
        while GOMIN:
            MIN_STAT = 0
            # MIN
            # TODO: bgoli: see whether this also works with 'minimize'
            glpk_setObjective(cpx, 'min%s' % R, [(1, R)], 'min', reset=True)
            ##  cplx_setBounds(c, id, min=None, max=None) # think about this
            MIN_STAT = glpk_Solve(cpx, method=method)
            if MIN_STAT == 'LPS_OPT':
                MIN_STAT = 1
            elif MIN_STAT == 'LPS_UNBND':
                MIN_STAT = 2
            elif MIN_STAT == 'LPS_NOFEAS':
                MIN_STAT = 3
            else:
                MIN_STAT = 4
            if MIN_STAT == 1: # solved
                min_oval = cpx.obj.value
            elif MIN_STAT == 2: # unbound
                min_oval = -numpy.Inf
            elif MIN_STAT == 3:
                #min_oval = pre_sol[R]
                min_oval = numpy.NaN
            else: # other failure
                min_oval = numpy.NaN
            if debug:
                cpx.write(cpxlp=os.path.join(debug_dir, '%smin.lp' % R))
            if MIN_STAT >= 3:
                if gomin_cntr == 0:
                    cpx.erase()
                    del cpx
                    time.sleep(0.1)
                    cpx = glpk.LPX(mps=mps_filename)
                    gomin_cntr += 1
                else:
                    GOMIN = False
            else:
                GOMIN = False
            if gomin_cntr >= max_error_iter:
                GOMIN = False

        GOMAX = True
        gomax_cntr = 0
        while GOMAX:
            MAX_STAT = 0
            # MAX
            glpk_setObjective(cpx, 'max%s' % R, expr=None, sense='max', reset=False)
            ##  cplx_setBounds(c, id, min=None, max=None) # think about this
            MAX_STAT = glpk_Solve(cpx, method=method)
            if MAX_STAT == 'LPS_OPT':
                MAX_STAT = 1
            elif MAX_STAT == 'LPS_UNBND':
                MAX_STAT = 2
            elif MAX_STAT == 'LPS_NOFEAS':
                MAX_STAT = 3
            else:
                MAX_STAT = 4
            if MAX_STAT == 1: # solved
                max_oval = cpx.obj.value
            elif MAX_STAT == 2: # unbound
                max_oval = numpy.Inf
            elif MAX_STAT == 3: # infeasable
                #max_oval = pre_sol[R]
                max_oval = numpy.NaN
            else: # other fail
                max_oval = numpy.NaN

            if MAX_STAT >= 3:
                if gomax_cntr == 0:
                    cpx.erase()
                    del cpx
                    time.sleep(0.1)
                    cpx = glpk.LPX(mps=mps_filename)
                    gomax_cntr += 1
                else:
                    GOMAX = False
            else:
                GOMAX = False
            if gomax_cntr >= max_error_iter:
                GOMAX = False

        # check for solver going berserk
        if MIN_STAT > 1 and MAX_STAT > 1:
            print(Ridx)
            time.sleep(1)

        # enables using the default value as a solution if the solver fails
        if pre_opt and default_on_fail:
            if MAX_STAT > 1 and not MIN_STAT > 1:
                max_oval = pre_sol[R]
            if MIN_STAT > 1 and not MAX_STAT > 1:
                min_oval = pre_sol[R]

        if debug:
            cpx.write(cpxlp=os.path.join(debug_dir, '%smax.lp' % R))

        OUTPUT_ARRAY[Ridx,0] = pre_sol[R]
        if R in REDUCED_COSTS:
            OUTPUT_ARRAY[Ridx,1] = REDUCED_COSTS[R]
        OUTPUT_ARRAY[Ridx,2] = min_oval
        OUTPUT_ARRAY[Ridx,3] = max_oval
        OUTPUT_ARRAY[Ridx,4] = round(abs(max_oval - min_oval), roundoff_span)
        OUTPUT_ARRAY[Ridx,5] = MIN_STAT
        OUTPUT_ARRAY[Ridx,6] = MAX_STAT
        if markupmodel:
            REAC = fba.getReaction(R)
            REAC.setValue(pre_sol[R])
            REAC.fva_min = min_oval
            REAC.fva_max = max_oval
            REAC.fva_status = (MIN_STAT, MAX_STAT)
            if R in REDUCED_COSTS:
                REAC.reduced_costs = REDUCED_COSTS[R]
        if not quiet and MAX_STAT > 1 or MIN_STAT > 1:
            print('Solver fail for reaction \"{}\" (MIN_STAT: {} MAX_STAT: {})'.format(R, MIN_STAT, MAX_STAT))
        cntr += 1
        if cntr == 200:
            tcnt += cntr
            print('FVA has processed {} of {} reactions'.format(tcnt, NUM_FLX))
            cntr = 0

    os.remove(mps_filename)

    #cpx.write(cpxlp='thefinaldebug.lp')
    del cpx
    print('\nSinglecore FVA took: {} min (1 process)\n'.format((time.time()-s2time)/60.))
    print('Output array has columns:')
    print('Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus')
    return OUTPUT_ARRAY, OUTPUT_NAMES


def glpk_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=False, oldlpgen=True, with_reduced_costs='unscaled', method='s'):
    """
    This is a utility function that does a presolve for FVA, MSAF etc. Generates properly formatted
    empty objects if pre_opt == False

     - *pre_opt* a boolean
     - *fba* a CBModel object
     - *objF2constr* add objective function as constraint
     - *quiet* [default=False] supress cplex output
     - *with_reduced_costs* [default='unscaled'] can be 'scaled' or 'unscaled'
     - *method* [default='s'] select the GLPK solver method, see the GLPK documentation for details

       - 's': simplex
       - 'i': interior
       - 'e': exact


    Returns: pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS

    """

    cpx = glpk_constructLPfromFBA(fba, fname=None)
    OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    if pre_opt:
        status = glpk_Solve(cpx, method=method)
        if glpk_getSolutionStatus(cpx) == 'LPS_OPT':
            print('Valid Presolution')
            OPTIMAL_PRESOLUTION = True
            pre_sol, pre_oid, pre_oval = glpk_getOptimalSolution(cpx)
            fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = pre_sol, pre_oval
            scaled = False
            if with_reduced_costs == 'scaled':
                REDUCED_COSTS = glpk_getReducedCosts(cpx, scaled=True)
            elif with_reduced_costs == 'unscaled':
                REDUCED_COSTS = glpk_getReducedCosts(cpx, scaled=False)
        else:
            print('Invalid Presolution, because {}'.format(glpk_getSolutionStatus(cpx) ))
            OPTIMAL_PRESOLUTION = False
            pre_sol = {}
            for r in fba.reactions:
                pre_sol.update({r.getPid : 0.0})
                r.reduced_cost = 0.0
            pre_oval = 0.0
            pre_oid = 'None'
            raise RuntimeError('\nPresolve failed to optimize this model and cannot continue!')
    else:
        pre_sol = {}
        for r in fba.reactions:
            pre_sol.update({r.getPid() : 0.0})
        if objF2constr:
            pre_oval = fba.objectives[fba.activeObjIdx].value
            pre_oid = fba.objectives[fba.activeObjIdx].getPid()
        else:
            pre_oval = 0.0
            pre_oid = 'None'
    for r in fba.reactions:
        r.value = pre_sol[r.id]
    return cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS

def glpk_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, oval, tol, optPercentage):
    """
    Take the objective function and "optimum" value and add it as a constraint
     - *cpx* a cplex object
     - *oval* the objective value
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)

    """

    # generate new constraint from old objective value (use non-zero coefficients)
    LCS = []
    LCN = []
    new_constraint = []
    for v_ in range(len(cpx.cols)):
        LCN.append(cpx.cols[v_].name)
        LCS.append(cpx.obj[v_])
        if LCS[-1] != 0.0:
            new_constraint.append((LCS[-1], LCN[-1]))

    if rhs_sense == 'equal':
        pass
    elif cpx.obj.maximize and (rhs_sense == 'upper'):
        print('\nWarning: RHS sense error: \"upper\" does not match \"maximize\" changing to \"lower\"')
        rhs_sense = 'lower'
        time.sleep(1)
    elif not cpx.obj.maximize and (rhs_sense == 'lower'):
        print('\nWarning: RHS sense error: \"lower\" does not match \"minimize\" changing to \"upper\"')
        rhs_sense = 'upper'
        time.sleep(1)
    else:
        print('\nRHS sense ok.')

    # set objective constraint
    if rhs_sense == 'equal':
        cplx_setSingleConstraint(cpx, 'ObjCstr', expr=new_constraint, sense='E', rhs=oval)
        ##  cplx_setSingleConstraint(cpx, 'ObjCstr', expr=[(1, pre_oid)], sense='E', rhs=oval)
    elif rhs_sense == 'upper':
        if tol != None:
            ub =  numpy.ceil(oval/tol)*tol*optPercentage/100.0
        else:
            ub = oval*(optPercentage/100.0)
        glpk_setSingleConstraint(cpx, 'ObjCstr', expr=new_constraint, sense='L', rhs=ub)
        ##  cplx_setSingleConstraint(cpx, 'ObjCstr', expr=[(1, pre_oid)], sense='L', rhs=ub)
    elif rhs_sense == 'lower':
        if tol != None:
            lb =  numpy.floor(oval/tol)*tol*optPercentage/100.0
        else:
            lb = oval*(optPercentage/100.0)
        glpk_setSingleConstraint(cpx, 'ObjCstr', expr=new_constraint, sense='G', rhs=lb)
        ##  cplx_setSingleConstraint(cpx, 'ObjCstr', expr=[(1, pre_oid)], sense='G', rhs=lb)
    else:
        raise RuntimeError("\nInvalid RHS sense: %s" % rhs_sense)

def glpk_setSingleConstraint(c, cid, expr=[], sense='E', rhs=0.0):
    """
    Sets a new sigle constraint to a GLPK model

     - *c* a GLPK instance
     - *cid* the constraint id
     - *expr* a list of (coefficient, name) pairs
     - *sense* [default='G'] LGE
     - *rhs* [default=0.0] the right hand side

    """

    baseRows = len(c.rows)
    c.rows.add(1)
    ofnames = [a[1] for a in expr]
    ofval = [a[0] for a in expr]

    # name and coefficients
    newCon = []
    for cidx_ in range(len(c.cols)):
        if c.cols[cidx_].name in ofnames:
            newCon.append((cidx_, ofval[ofnames.index(c.cols[cidx_].name)]))
            print(newCon[-1])
        else:
            newCon.append((cidx_, 0.0))
    c.rows[baseRows].name = cid
    c.rows[baseRows].matrix = newCon

    # sense and rhs
    if sense in ['<=','<','L']:
        c.rows[baseRows].bounds = None, rhs
    elif sense in ['>=','>','G']:
        c.rows[baseRows].bounds = rhs, None
    elif sense in ['=','E']:
        c.rows[baseRows].bounds = rhs
    else:
        raise RuntimeError('INFO: invalid operator: %s' % sense)
    #c.write(cpxlp='test.setc.lp')

def glpk_setObjective(c, oid, expr=None, sense='maximize', reset=True):
    """
    Set a new objective function note that there is a major memory leak in
    `c.variables.get_names()` whch is used when reset=True. If this is a problem
    use cplx_setObjective2 which takes *names* as an input:

     - *c* a GLPK LP object
     - *oid* the r_id of the flux to be optimized
     - *expr* a list of (coefficient, flux) pairs
     - *sense* 'maximize'/'minimize'
     - *reset* [default=True] reset all objective function coefficients to zero

    """
    sense = sense.lower()
    if sense == 'max': sense = 'maximize'
    if sense == 'min': sense = 'minimize'
    if sense in ['maximise', 'minimise']:
        sense = sense.replace('se','ze')
    assert sense in ['maximize', 'minimize'], "\nsense must be ['maximize', 'minimize'] not %s" % sense

    if expr != None:
        ofnames = [a[1] for a in expr]
        ofval = [a[0] for a in expr]

        for cidx_ in range(len(c.cols)):
            if c.cols[cidx_].name in ofnames:
                c.obj[cidx_] = ofval[ofnames.index(c.cols[cidx_].name)]
            elif reset:
                c.obj[cidx_] = 0.0
            else:
                c.obj[cidx_] = c.obj[cidx_]
    c.obj.name = oid

    if sense == 'minimize':
        c.obj.maximize = False
        if __DEBUG__: print('Set minimizing')
    else:
        c.obj.maximize = True
        if __DEBUG__: print('Set maximizing')
    #c.write(cpxlp='test_obja.lp')

def getReducedCosts(fba):
    """
    Get a dictionary of reduced costs for each reaction/flux

    """
    output = {}
    for r in fba.reactions:
        output.update({r.getPid() : r.reduced_cost})
    return output

def setReducedCosts(fba, reduced_costs):
    """
    For each reaction/flux, sets the attribute "reduced_cost" from a dictionary of
    reduced costs

     - *fba* an fba object
     - *reduced_costs* a dictionary of {reaction : value} pairs

    """
    if len(reduced_costs) == 0:
        pass
    else:
        for r in fba.reactions:
            if r.getPid() in reduced_costs:
                r.reduced_cost = reduced_costs[r.getPid()]
            else:
                r.reduced_cost = None


#def cplx_getDualValues(c):
    #"""
    #Get the get the dual values of the solution

     #- *c* a CPLEX LP

    #Output is a dictionary of {name : value} pairs

    #"""
    #d_names = c.linear_constraints.get_names()
    #d_values = c.solution.get_dual_values()
    #output = {}
    #for j in range(len(d_names)):
        #output.update({d_names[j] : d_values[j]})
    #return output

#def cplx_getSensitivities(c):
    #"""
    #Get the sensitivities of each constraint on the objective function with input

     #- *c* a CPLEX LP

    #Output is a tuple of bound and objective sensitivities where the objective
    #sensitivity is described in the CPLEX reference manual as::

     #... the objective sensitivity shows each variable, its reduced cost and the range over
     #which its objective function coefficient can vary without forcing a change
     #in the optimal basis. The current value of each objective coefficient is
     #also displayed for reference.

     #- *objective coefficient sensitivity* {flux : (reduced_cost, lower_obj_sensitivity, coeff_value, upper_obj_sensitivity)}
     #- *rhs sensitivity* {constraint : (low, value, high)}
     #- *bound sensitivity ranges* {flux : (lb_low, lb_high, ub_low, ub_high)}

    #"""
    #SENSE_RHS = {}
    #SENSE_BND = {}
    #SENSE_OBJ = {}
    #c_names = c.linear_constraints.get_names()
    #rhs_val = c.linear_constraints.get_rhs()
    #j_names = c.variables.get_names()

    #rhs_sense = c.solution.sensitivity.rhs()
    #bnd_sense = c.solution.sensitivity.bounds()
    #obj_sense = c.solution.sensitivity.objective()
    #obj_coeff = c.objective.get_linear()
    #red_cost = c.solution.get_reduced_costs()

    #for r in range(c.variables.get_num()):
        #SENSE_BND.update({j_names[r] : (bnd_sense[r][0], bnd_sense[r][1], bnd_sense[r][2], bnd_sense[r][3])})
        #SENSE_OBJ.update({j_names[r] : (red_cost[r], obj_sense[r][0], obj_coeff[r], obj_sense[r][1])})

    #for s in range(c.linear_constraints.get_num()):
        #SENSE_RHS.update({c_names[s] : (rhs_sense[s][0], rhs_val[s], rhs_sense[s][1])})

    #return (SENSE_OBJ, SENSE_RHS, SENSE_BND)


def glpk_MinimizeSumOfAbsFluxes(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=False, debug=False, objective_coefficients={}, return_lp_obj=False, oldlpgen=False, with_reduced_costs=None, method='s'):
    """
    Minimize the sum of absolute fluxes sum(abs(J1) + abs(J2) + abs(J3) ... abs(Jn)) by adding two constraints per flux
    and a variable representing the absolute value:

     Min: Ci abs_Ji
      Ji - abs_Ji <= 0
      Ji + abs_Ji >= 0

     Such that:
      NJi = 0
      Jopt = opt

    returns the value of the flux minimization objective function (not the model objective function which remains unchanged from)

    Arguments:

     - *fba* an FBA model object
     - *selected reactions* [default=None] means use all reactions otherwise use the reactions listed here
     - *pre_opt* [default=True] attempt to presolve the FBA and report its results in the ouput, if this is disabled and *objF2constr* is True then the vid/value of the current active objective is used
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)
     - *work_dir* [default=None] the MSAF working directory for temporary files default = cwd+fva
     - *debug* [default=False] if True write out all the intermediate MSAF LP's into work_dir
     - *quiet* [default=False] if enabled supress CPLEX output
     - *objF2constr* [default=True] add the model objective function as a constraint using rhs_sense etc. If
       this is True with pre_opt=False then the id/value of the active objective is used to form the constraint
     - *objective_coefficients* [default={}] a dictionary of (reaction_id : float) pairs that provide the are introduced as objective coefficients to the absolute flux value. Note that the default value of the coefficient (non-specified) is +1.
     - *return_lp_obj* [default=False] off by default when enabled it returns the CPLEX LP object
     - *with_reduced_costs* [default=None] if not None should be 'scaled' or 'unscaled'
     - *method* [default='s'] select the GLPK solver method, see the GLPK documentation for details

       - 's': simplex
       - 'i': interior
       - 'e': exact


    With outputs:

     - *fba* an update instance of a CBModel. Note that the FBA model objective function value is the original value set as a constraint

    """

    if with_reduced_costs == 'scaled':
        fba.SCALED_REDUCED_COSTS = True
    elif with_reduced_costs == 'unscaled':
        fba.SCALED_REDUCED_COSTS = False
    elif fba.SCALED_REDUCED_COSTS == True:
        with_reduced_costs = 'scaled'
    elif fba.SCALED_REDUCED_COSTS == False:
        with_reduced_costs = 'unscaled'
    else:
        with_reduced_costs = None

    print('\nINFO: using \"{}\" reduced costs.\n'.format(with_reduced_costs))

    if work_dir == None:
        work_dir = os.getcwd()
    else:
        assert os.path.exists(work_dir), '\nWhat did you think would happen now!'
    if debug:
        debug_dir = os.path.join(work_dir,'DEBUG')
        if not os.path.exists(debug_dir):
            os.mkdir(debug_dir)

    # generate a presolution
    print('GLPK is using solver option: "{}"'.format(method))
    cpx = OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = glpk_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, with_reduced_costs=with_reduced_costs, method=method)
    # if required add the objective function as a constraint
    if objF2constr:
        glpk_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, tol, optPercentage)

    STORED_OPT = None
    if pre_opt:
        STORED_OPT = pre_oval
    else:
        STORED_OPT = fba.getActiveObjective().getValue()

    # minimize the absolute sum of fluxes
    if selected_reactions != None:
        rids = fba.getReactionIds()
        for r in selected_reactions:
            assert r in rids, "\n%s is not a valid reaction name" % r
    else:
        selected_reactions = fba.getReactionIds()

    # this removes the objective function
    ##  fba_obj_ids = fba.getActiveObjective().getFluxObjectiveReactions()
    ##  print fba_obj_ids
    ##  for R in fba_obj_ids:
        ##  if R in selected_reactions:
            ##  selected_reactions.pop(selected_reactions.index(R))
    ##  del R, fba_obj_ids
    NUM_FLX = len(selected_reactions)
    print('Total number of reactions: {}'.format(NUM_FLX))
    abs_selected_reactions = ['abs_%s' % r for r in selected_reactions]
    # absJ implicitly assumed to be absJ >= 0
    #C#cpx.variables.add(names=abs_selected_reactions)
    basecols = len(cpx.cols)
    cpx.cols.add(len(abs_selected_reactions))
    for v_ in range(len(abs_selected_reactions)):
        cpx.cols[basecols+v_].name = abs_selected_reactions[v_]
        cpx.cols[basecols+v_].bounds = 0, None

    """
     J - abs_J <= 0
     J + abs_J >= 0
    """
    baserows = len(cpx.rows)
    #tnew = time.time()
    rowNames = tuple([c.name for c in cpx.rows])
    colNames = tuple([c.name for c in cpx.cols])
    cRange = range(baserows, baserows+(len(selected_reactions)))
    rcntr = 0
    cpx.rows.add(len(selected_reactions)*2)
    for r_ in cRange:
        # name and coefficients
        cpx.rows[r_].name = 'abs{}a'.format(selected_reactions[rcntr])
        cpx.rows[r_].matrix = [(colNames.index(selected_reactions[rcntr]), 1),\
                               (colNames.index(abs_selected_reactions[rcntr]), -1)]
        cpx.rows[r_].bounds = None, 0
        cpx.rows[r_+len(selected_reactions)].name = 'abs{}b'.format(selected_reactions[rcntr])
        cpx.rows[r_+len(selected_reactions)].matrix = [(colNames.index(selected_reactions[rcntr]), 1),\
                                (colNames.index(abs_selected_reactions[rcntr]), 1)]
        cpx.rows[r_+len(selected_reactions)].bounds = 0, None
        rcntr += 1

    glpk_setObjective(cpx, 'MAFS', [(1, j) for j in abs_selected_reactions], 'min', reset=True)

    if debug:
        cpx.write(cpxlp=os.path.join(debug_dir, 'MSAF_base_(%s).lp' % time.time()))
        #cplx_writeLPtoLPTfile(cpx, , title=None, Dir=debug_dir)
        ##  cplx_writeLPtoLPTfile(cpx, 'MSAF_base_%s' % time.time() , title=None, Dir=debug_dir)

    glpk_Solve(cpx, method=method)
    glpk_setFBAsolutionToModel(fba, cpx, with_reduced_costs=with_reduced_costs)

    #cpx.write(cpxlp='test.lp'); return cpx

    minSum = cpx.obj.value
    fba.setAnnotation('min_flux_sum', minSum)
    fba.getActiveObjective().setValue(STORED_OPT)
    print('\nMinimizeSumOfAbsFluxes objective value: {}\n'.format(minSum))
    if quiet:
        pass
        #cplx_setOutputStreams(cpx, mode='default')
    if not return_lp_obj:
        return minSum
    else:
        return cpx

def glpk_writeLPtoLPTfile(c, filename, title=None, Dir=None):
    """
    Write out a GLPK model as an LP format file

    """
    if Dir != None:
        filename = os.path.join(Dir, filename)
    #if title != None:
        #c.set_problem_name(title)
    #c.write(filename+'.lp', filetype='lp')
    c.write(cpxlp=filename+'.lp')
    print('LP output as {}'.format(filename+'.lp'))
