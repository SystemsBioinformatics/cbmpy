"""
CBMPy: CBCPLEX module
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2016 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBCPLEX.py 569 2017-03-10 15:50:07Z bgoli $)

"""
## gets rid of "invalid variable name" info
# pylint: disable=C0103
## gets rid of "line to long" info
# pylint: disable=C0301
## use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, gc
# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio

HAVE_SYMPY = False
try:
    import sympy
    if int(sympy.__version__.split('.')[1]) >= 7 and int(sympy.__version__.split('.')[2]) >= 4:
        HAVE_SYMPY = True
    else:
        del sympy
except ImportError:
    HAVE_SYMPY = False
HAVE_SCIPY = False
try:
    from scipy.sparse import csr
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

import numpy
import cplex
#from cplex.exceptions import CplexError
from . import CBWrite, CBTools
from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


##  class FluxSensitivity(object):
    ##  """
    ##  Contains variable (Flux) sensitivities and reduced cost
    ##  """
    ##  flux = None
    ##  lower_bound_low = None
    ##  lower_bound_high = None
    ##  lower_bound_value = None
    ##  upper_bound_low = None
    ##  upper_bound_high = None
    ##  upper_bound_value = None
    ##  objective_coefficient_value = None
    ##  objective_coefficient_low = None
    ##  objective_coefficient_high = None
    ##  reduced_cost = None

    ##  def getUpperBoundRange(self):
        ##  """
        ##  Returns the absolute range of the upper bound.
        ##  """
        ##  return abs(self.upper_bound_high-self.upper_bound_low)

    ##  def getLowerBoundRange(self):
        ##  """
        ##  Returns the absolute range of the lower bound.
        ##  """
        ##  return abs(self.lower_bound_high-self.lower_bound_low)

    ##  def getObjectiveCoefficientRange(self):
        ##  """
        ##  Returns the absolute range of the objective flux coefficient.
        ##  """
        ##  return abs(self.objective_coefficient_high-objective_coefficient_low)

##  class RHSsensitivity(object):
    ##  """
    ##  Class containing RHS sensitivities
    ##  """
    ##  species = None
    ##  rhs_low = None
    ##  rhs_high = None
    ##  rhs_value = None

    ##  def getRHSRange(self):
        ##  """
        ##  Returns the absolute range of the objective flux coefficient.
        ##  """
        ##  return abs(self.rhs_high-rhs_low)

##  class Sensitivities(object):
    ##  """
    ##  Container class that holds sensitivities
    ##  """
    ##  flux_ids = None
    ##  species_ids = None

    ##  def __init__(self):
        ##  self.flux_ids = []
        ##  self.species_ids = []

    ##  def addFluxSensitivity(self, fsense):
        ##  self.__setattr__(fsense.flux, fsense)
        ##  self.flux_ids.append(fsense.flux)

    ##  def addRHSsensitivity(self, rhsense):
        ##  self.__setattr__(rhsense.species, rhsense)
        ##  self.species_ids.append(rhsense.species)


"""
Sets the noise level of the solver CPLX_RESULT_STREAM can be::

 - *None* silent i.e. no output
 - *file* set solver to silent and output logs to *CPLX_RESULT_STREAM_FILE* cplex_output.log
 - *iostream* set solver to silent and output logs to *CPLX_RESULT_STREAM_IO* csio
 - *'default'* or anything else noisy with full output closes STREAM_IO and STREAM_FILE (default)
"""
#CPLX_RESULT_STREAM, CPLX_RESULT_STREAM_FILE, CPLX_RESULT_STREAM_IO, CPLX_SILENT_MODE
CPLX_RESULT_STREAM = 'default'
CPLX_RESULT_STREAM_FILE = None
CPLX_RESULT_STREAM_IO = None
CPLX_SILENT_MODE = True

CPLX_LP_PARAMETERS = {
    'simplex.tolerances.optimality' : 1e-6,
    'simplex.tolerances.feasibility' : 1e-6
    }


def cplx_fixConSense(operator):
    """
    Fixes the sense of inequality operators, returns corrected sense symbol

    - *operator* the operator to check

    """

    if operator in ['<=', '<','L']:
        oPr = 'L'
    elif operator in ['>=', '>','G']:
        oPr = 'G'
    elif operator in ['=', 'E']:
        oPr = 'E'
    else:
        raise RuntimeError('INFO: invalid operator: %s' % operator)
    return oPr


def cplx_constructLPfromFBA(fba, fname=None):
    '''
    Create a CPLEX LP in memory.
    - *fba* an FBA object
    - *fname* optional filename if defined writes out the constructed lp

    '''
    _Stime = time.time()
    # defines
    #csense = 'E' # constraint sense

    # define model and add variables
    lp = cplex.Cplex()
    # define simplex tolerances for the model
    lp.parameters.simplex.tolerances.optimality.set(CPLX_LP_PARAMETERS['simplex.tolerances.optimality'])
    lp.parameters.simplex.tolerances.feasibility.set(CPLX_LP_PARAMETERS['simplex.tolerances.feasibility'])
    #print(lp.parameters.simplex.get_changed())
    lp.set_problem_name('%s' % (fba.getId()))
    lp.variables.add(names=fba.N.col)
    try:
        # define objective
        osense = fba.getActiveObjective().operation.lower()
        if osense in ['minimize', 'minimise', 'min']:
            lp.objective.set_sense(lp.objective.sense.minimize)
        elif osense in ['maximize', 'maximise', 'max']:
            lp.objective.set_sense(lp.objective.sense.maximize)
        else:
            raise RuntimeError('\n%s - is not a valid objective operation' % osense)
        lp.objective.set_name(fba.getActiveObjective().getId())
        lp.objective.set_linear([(fo.reaction, fo.coefficient) for fo in fba.getActiveObjective().fluxObjectives])
    except AttributeError:
        print('\nWARNING(CPLEX create LP): no objective function defined')

    ## old fashioned way
    #told = time.time()
    #lin_expr = []
    #rhs =[]
    #names = []
    #senses = []
    #for lc_ in range(fba.N.shape[0]):
        #ids = []
        #val = []
        #for c_ in range(len(fba.N.array[lc_])):
            #if fba.N.array[lc_, c_] != 0.0:
                #ids.append(fba.N.col[c_])
                #val.append(fba.N.array[lc_, c_])
        #lin_expr.append(cplex.SparsePair(ids, val))
        #rhs.append(fba.N.RHS[lc_])
        #senses.append(fba.N.operators[lc_])
        #names.append(fba.N.row[lc_])
    #lp.linear_constraints.add(lin_expr=lin_expr, senses=senses, rhs=rhs, names=names)
    #print 'Old style lc:', time.time() - told

    ## the numpy way
    #tnew = time.time()
    lin_expr = []
    rhs = []
    names = []
    senses = []
    if HAVE_SYMPY and fba.N.__array_type__ == sympy.MutableDenseMatrix:
        print('INFO: CPLEX requires floating point, converting N')
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

    GOSCI = False
    if HAVE_SCIPY and fba.N.__array_type__ == csr.csr_matrix:
        GOSCI = True

    for n in range(Nmat.shape[0]):
        if not GOSCI:
            nz = Nmat[n].nonzero()[0]
        else:
            nz = Nmat[n].nonzero()[1]
        lin_expr.append(cplex.SparsePair([fba.N.col[c] for c in nz], [Nmat[n, c] for c in nz]))
        rhs.append(RHSmat[n])
        senses.append(cplx_fixConSense(fba.N.operators[n]))
        names.append(fba.N.row[n])
    #print senses
    lp.linear_constraints.add(lin_expr=lin_expr, senses=senses, rhs=rhs, names=names)
    #print 'New style lc:', time.time() - tnew

    # add user defined constraints
    if fba.CM != None:
        ## the numpy way
        #t2new = time.time()
        lin_expr = []
        rhs = []
        names = []
        senses = []
        for n in range(CMmat.shape[0]):
            if not GOSCI:
                nz = CMmat[n].nonzero()[0]
            else:
                nz = CMmat[n].nonzero()[1]
            lin_expr.append(cplex.SparsePair([fba.CM.col[c] for c in nz], [CMmat[n, c] for c in nz]))
            rhs.append(CMrhs[n])
            senses.append(cplx_fixConSense(fba.CM.operators[n]))
            names.append(fba.CM.row[n])
        lp.linear_constraints.add(lin_expr=lin_expr, senses=senses, rhs=rhs, names=names)
        #print 'New style lc:', time.time() - t2new

    # add bounds
    lb = []
    ub = []
    for b_ in fba.flux_bounds:
        btype = b_.getType()
        bvalue = b_.getValue()
        if bvalue in ['Infinity', 'inf', 'Inf', 'infinity']:
            bvalue = cplex.infinity
        elif bvalue in ['-Infinity', '-inf', '-Inf', '-infinity']:
            bvalue = -cplex.infinity
        elif numpy.isinf(bvalue):
            if bvalue > 0.0:
                bvalue = cplex.infinity
            else:
                bvalue = -cplex.infinity
        if btype == 'lower':
            lb.append((b_.reaction, bvalue))
        elif btype == 'upper':
            ub.append((b_.reaction, bvalue))
        elif btype == 'equality':
            ub.append((b_.reaction, bvalue))
            lb.append((b_.reaction, bvalue))
    if len(lb) > 0:
        lp.variables.set_lower_bounds(lb)
    if len(ub) > 0:
        lp.variables.set_upper_bounds(ub)

    print('\ncplx_constructLPfromFBA time: {}\n'.format(time.time() - _Stime))
    if fname != None:
        lp.write(fname+'.lp', filetype='lp')
    return lp

def cplx_analyzeModel(f, lpFname=None, return_lp_obj=False, with_reduced_costs='unscaled', with_sensitivity=False,\
                      del_intermediate=False, build_n=True, quiet=False, oldlpgen=False, method='o'):
    """
    Optimize a model and add the result of the optimization to the model object
    (e.g. `reaction.value`, `objectiveFunction.value`). The stoichiometric
    matrix is automatically generated. This is a common function available
    in all solver interfaces. By default returns the objective function value

     - *f* an instantiated PySCeSCBM model object
     - *lpFname* [default=None] the name of the intermediate LP file. If not specified no LP file is produced
     - *return_lp_obj* [default=False] off by default when enabled it returns the CPLEX LP object
     - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
       or anything else which is interpreted as 'None'. Scaled means s_rcost = (r.reduced_cost*rval)/obj_value
     - *with_sensitivity* [default=False] add solution sensitivity information (not yet implemented)
     - *del_intermediate* [default=False] redundant except if output file is produced and deleted (not useful)
     - *build_n* [default=True] generate stoichiometry from the reaction network (reactions/reagents/species)
     - *quiet* [default=False] suppress cplex output
     - *method* [default='o'] choose the CPLEX method to use for solution, default is automatic. See CPLEX reference manual for details

       - 'o': auto
       - 'p': primal
       - 'd': dual
       - 'b': barrier (no crossover)
       - 'h': barrier
       - 's': sifting
       - 'c': concurrent


    """

    if build_n:
        f.buildStoichMatrix()

    fid = f.id

    if with_reduced_costs == 'scaled':
        f.SCALED_REDUCED_COSTS = True
    elif with_reduced_costs == 'unscaled':
        f.SCALED_REDUCED_COSTS = False
    else:
        f.SCALED_REDUCED_COSTS = None

    if lpFname == None:
        ##  lpFname = '_cplxtmp_%s_.tmp' % random.randint(0,100)
        lpFname = '_cplxtmp_.tmp'
    f.id = lpFname
    _Stime = time.time()
    if oldlpgen:
        LPF = CBWrite.writeModelLP(f, quiet=quiet)
        LPF = os.path.abspath(LPF)
        flp = cplx_getModelFromLP(LPF)
    else:
        if lpFname == '_cplxtmp_.tmp':
            flp = cplx_constructLPfromFBA(f, fname=None)
        else:
            flp = cplx_constructLPfromFBA(f, fname=lpFname)
    if not quiet:
        print('\ncplx_analyzeModel FBA --> LP time: {}\n'.format(time.time() - _Stime))
    f.id = fid
    if not quiet and CPLX_SILENT_MODE:
        cplx_setOutputStreams(flp, mode='default')
    elif quiet and CPLX_SILENT_MODE:
        pass
    elif quiet:
        cplx_setOutputStreams(flp, mode=None)
    cplx_Solve(flp, method)
    cplx_setFBAsolutionToModel(f, flp, with_reduced_costs=with_reduced_costs)
    cplx_setSolutionStatusToModel(f, flp)
    if oldlpgen and del_intermediate:
        os.remove(LPF)
    objv = f.getActiveObjective().getValue()
    if not quiet:
        print('\nanalyzeModel objective value: {}\n'.format(objv))
    if return_lp_obj:
        return flp
    else:
        del flp
        return objv

def cplx_setFBAsolutionToModel(fba, lp, with_reduced_costs='unscaled'):
    """
    Sets the FBA solution from a CPLEX solution to an FBA object

     - *fba* and fba object
     - *lp* a CPLEX LP object
     - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
       or anything else which is interpreted as None. Scaled is: s_rcost = (r.reduced_cost*rval)/obj_value

    """
    sol, objname, objval = cplx_getOptimalSolution(lp)
    if not CPLX_SILENT_MODE:
        print('Status:', cplx_getSolutionStatus(lp))
    if cplx_getSolutionStatus(lp) == 'LPS_OPT':
        fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = sol, objval
    else:
        fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = sol, numpy.NaN
    for r in fba.reactions:
        rid = r.getId()
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
        RC = cplx_getReducedCosts(lp, scaled=scaled)
        setReducedCosts(fba, RC)
    else:
        setReducedCosts(fba, {})

def getReducedCosts(fba):
    """
    Get a dictionary of reduced costs for each reaction/flux

    """
    output = {}
    for r in fba.reactions:
        output.update({r.getId() : r.reduced_cost})
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
            if r.getId() in reduced_costs:
                r.reduced_cost = reduced_costs[r.getId()]
            else:
                r.reduced_cost = None

def cplx_setOutputStreams(lp, mode='default'):
    """
    Sets the noise level of the solver, mode can be one of:

     - *None* silent i.e. no output
     - *'file'* set solver to silent and output logs to *CPLX_RESULT_STREAM_FILE* cplex_output.log
     - *'iostream'* set solver to silent and output logs to *CPLX_RESULT_STREAM_IO* csio
     - *'default'* or anything else noisy with full output closes STREAM_IO and STREAM_FILE (default)
    """
    global CPLX_RESULT_STREAM
    global CPLX_SILENT_MODE
    global CPLX_RESULT_STREAM_FILE
    global CPLX_RESULT_STREAM_IO
    CPLX_RESULT_STREAM = mode

    if CPLX_RESULT_STREAM == None:
        lp.set_log_stream(None)
        lp.set_results_stream(None)
        lp.set_warning_stream(None)
        CPLX_SILENT_MODE = True
    elif CPLX_RESULT_STREAM == 'file':
        if CPLX_RESULT_STREAM_FILE == None:
            CPLX_RESULT_STREAM_FILE = file('cplex_output.log', 'w')
        lp.set_log_stream(CPLX_RESULT_STREAM_FILE)
        lp.set_results_stream(CPLX_RESULT_STREAM_FILE)
        CPLX_SILENT_MODE = True
    elif CPLX_RESULT_STREAM == 'iostream':
        if CPLX_RESULT_STREAM_IO == None:
            CPLX_RESULT_STREAM_IO = csio.StringIO()
        lp.set_log_stream(CPLX_RESULT_STREAM_IO)
        lp.set_results_stream(CPLX_RESULT_STREAM_IO)
        CPLX_SILENT_MODE = True
    else:
        CPLX_RESULT_STREAM = 'default'
        CPLX_SILENT_MODE = False
        try:
            CPLX_RESULT_STREAM_IO.close()
        except AttributeError:
            pass
        try:
            CPLX_RESULT_STREAM_FILE.close()
        except AttributeError:
            pass

        CPLX_RESULT_STREAM_IO = None
        CPLX_RESULT_STREAM_FILE = None
        lp.set_log_stream(os.sys.stdout)
        lp.set_results_stream(os.sys.stdout)

def cplx_getModelFromLP(lptFile, Dir=None):
    """
    Load a LPT (CPLEX format) file and return a CPLX LP model

     - *lptfile* an CPLEX LP format file
     - *Dir* an optional directory

    """
    if Dir != None:
        assert os.path.exists(Dir), '\nIncorrect path'
        lptFile = os.path.join(Dir, lptFile)
    lp = cplex.Cplex(lptFile)
    # define simplex tolerances for the model
    lp.parameters.simplex.tolerances.optimality.set(CPLX_LP_PARAMETERS['simplex.tolerances.optimality'])
    lp.parameters.simplex.tolerances.feasibility.set(CPLX_LP_PARAMETERS['simplex.tolerances.feasibility'])
    return lp

def cplx_getCPLEXModelFromLP(lptFile, Dir=None):
    """
    Load a LPT (CPLEX format) file and return a CPLX LP model

     - *lptfile* an CPLEX LP format file
     - *Dir* an optional directory

    """
    return cplx_getModelFromLP(lptFile, Dir)


def cplx_getModelFromObj(fba):
    """
    Return a CPLEX object from a FBA model object (via LP file)
    """
    _Stime = time.time()
    modId = fba.id
    ##  fba.id = '_cplxtmp_%s_.tmp' % random.randint(0,5)
    fba.id = '_cplxtmp_.tmp'
    LPF = CBWrite.writeModelLP(fba, work_dir=os.getcwd())
    fba.id = modId
    x = cplx_getModelFromLP(LPF)
    print('\ncplx_getModelFromObj time: {}\n'.format(time.time() - _Stime))
    return x

def cplx_Solve(c, method='o'):
    alg = c.parameters.lpmethod.values

    if method == "o":
        c.parameters.lpmethod.set(alg.auto)
    elif method == "p":
        c.parameters.lpmethod.set(alg.primal)
    elif method == "d":
        c.parameters.lpmethod.set(alg.dual)
    elif method == "b":
        c.parameters.lpmethod.set(alg.barrier)
        c.parameters.barrier.crossover.set(c.parameters.barrier.crossover.values.none)
    elif method == "h":
        c.parameters.lpmethod.set(alg.barrier)
    elif method == "s":
        c.parameters.lpmethod.set(alg.sifting)
    elif method == "c":
        c.parameters.lpmethod.set(alg.concurrent)
    else:
        print("Unrecognized option, using automatic")
        c.parameters.lpmethod.set(alg.auto)
    try:
        c.solve()
    except cplex.exceptions.CplexSolverError as ex:
        print("Exception raised during solve:\n\t\"{}\"".format(ex))
        return None
    # solution.get_status() returns an integer code
    status = c.solution.get_status()
    if status == c.solution.status.optimal:
        if not CPLX_SILENT_MODE:
            print("INFO: Model is optimal:", status)
    elif status == c.solution.status.unbounded:
        if not CPLX_SILENT_MODE:
            print("INFO: Model is unbounded:", status)
        return status
    elif status == c.solution.status.infeasible:
        if not CPLX_SILENT_MODE:
            print("INFO: Model is infeasible:", status)
        return status
    elif status == c.solution.status.infeasible_or_unbounded:
        if not CPLX_SILENT_MODE:
            print("INFO: Model is infeasible or unbounded:", status)
        return status

    if not CPLX_SILENT_MODE:
        s_method = c.solution.get_method()
        print("Solution status = ", status, ":", end=" ")
        # the following line prints the status as a string
        print(c.solution.status[status])
        print("Solution method = ", s_method, ":", end=" ")
        print(c.solution.method[s_method])

    s_type = c.solution.get_solution_type()
    if s_type == c.solution.type.none:
        print("INFO: No solution available")
        return None
    else:
        if not CPLX_SILENT_MODE:
            print("Objective value = ", c.solution.get_objective_value())
        return status

def cplx_setMIPGapTolerance(c, tol):
    """
    Sets the the relative MIP gap tolerance
    """
    c.parameters.mip.tolerances.mipgap.set(tol)
    print('New MIP gap tolerance', c.parameters.mip.tolerances.mipgap.get())

# TODO
#def cplx_setSolverParameters(c, parset):
    #pass

def cplx_SolveMILP(c, auto_mipgap=False):
    """
    Solve and MILP

     - *auto_mipgap* auto decrease mipgap until mipgap == absmipgap

    """
    try:
        c.solve()
        if c.solution.get_status() == c.solution.status.optimal_tolerance:
            if auto_mipgap:
                while c.parameters.mip.tolerances.mipgap.get() >= c.parameters.mip.tolerances.absmipgap.get():
                    c.parameters.mip.tolerances.mipgap.set(c.parameters.mip.tolerances.mipgap.get()/10.0)
                    c.solve()
                    if c.solution.get_status() == c.solution.status.MIP_optimal:
                        break
                print('\n\nMILP solution gap tolerance set to: %s\n\n' % c.parameters.mip.tolerances.mipgap.get())
                time.sleep(5)
            if c.solution.get_status() == c.solution.status.optimal_tolerance:
                print('\n\nMILP solution found with gap tolerance consider adjusting gap tolerance with cplx_setGapTolerance\n\n')
                #raise RuntimeWarning
    except cplex.exceptions.CplexSolverError as ex:
        print("\nException raised during solve\n")
        print(ex)

def cplx_getOptimalSolution(c):
    """
    From a CPLX model extract a tuple of solution, ObjFuncName and ObjFuncVal
    """
    s_val = []
    s_name = []
    fba_sol = {}
    objf_name = None
    objf_val = None
    try:
        objf_name = c.objective.get_name()
        objf_val = c.solution.get_objective_value()
        s_val = c.solution.get_values()
        s_name = c.variables.get_names()
        for n in range(len(s_name)):
            fba_sol.update({s_name[n] : s_val[n]})
    except Exception as ex:
        print(ex)
        print('WARNING: No solution to get')
        s_val = []
        s_name = []
        fba_sol = {}
        objf_val = None
    del s_name, s_val
    return fba_sol, objf_name, objf_val

def cplx_getOptimalSolution2(c, names):
    """
    From a CPLX model extract a tuple of solution, ObjFuncName and ObjFuncVal
    """
    s_val = []
    s_name = []
    fba_sol = None
    objf_name = None
    objf_val = None
    try:
        s_val = c.solution.get_values()
        objf_name = c.objective.get_name()
        objf_val = c.solution.get_objective_value()
        s_name = names
    except AttributeError:
        print('No solution attribute on CPLEX object')
    fba_sol = {}
    for n in range(len(s_name)):
        fba_sol.update({s_name[n] : s_val[n]})
    del names, s_val
    return fba_sol, objf_name, objf_val

def cplx_writeLPsolution(fba_sol, objf_name, fname, Dir=None, separator=','):
    """
    This function writes the optimal solution, produced wth `cplx_getOptimalSolution` to file

     - *fba_sol* a dictionary of Flux : value pairs
     - *objf_name* the objective function flux id
     - *fname* the output filename
     - *Dir* [default=None] use directory if not None
     - *separator* [default=','] the column separator

    """
    if Dir != None:
        assert os.path.exists(Dir), '\nPath does not exist'
        fname = os.path.join(Dir, fname)
    fname += '.csv'
    F = file(fname, 'w')
    F.write('%s%s%s\n' % ('ObjectiveFunction', separator, objf_name))
    for r in fba_sol:
        F.write('%s%s%f\n' % (r, separator, fba_sol[r]))
    F.flush()
    F.close()
    print('CSV exported to {}'.format(fname))

def cplx_getShadowPrices(c):
    """
    Returns a dictionary of shadow prices containing 'Rid' : (lb, rhs, ub)

     - *c* a cplex LP object

    """
    c_names = c.linear_constraints.get_names()
    rhs_sense = c.solution.sensitivity.rhs()
    rhs_val = c.linear_constraints.get_rhs()
    output = {}
    for s in range(c.linear_constraints.get_num()):
        output.update({c_names[s] : (rhs_sense[s][0], rhs_val[s], rhs_sense[s][1])})
    return output

def cplx_getReducedCosts(c, scaled=False):
    """
    Extract ReducedCosts from LP and return as a dictionary 'Rid' : reduced cost

     - *c* a cplex LP object
     - *scaled* scale the reduced cost by the optimal flux value

    """
    s_name = c.variables.get_names()
    r_costs = c.solution.get_reduced_costs()
    objf_val = c.solution.get_objective_value()
    output = {}
    s_val = None
    if scaled:
        s_val = c.solution.get_values()
    for v in range(len(s_name)):
        if scaled:
            try:
                r_val = r_costs[v]*s_val[v]/objf_val
            except Exception:
                r_val = 0.0
        else:
            r_val = r_costs[v]
        output.update({s_name[v] : r_val})
    del s_name, r_costs, s_val
    return output

def cplx_getSensitivities(c):
    """
    Get the sensitivities of each constraint on the objective function with inpt

     - *c* a CPLEX LP

    Output is a tuple of bound and objective sensitivities where the objective
    sensitivity is described in the CPLEX reference manual as::

     ... the objective sensitivity shows each variable, its reduced cost and the range over
     which its objective function coefficient can vary without forcing a change
     in the optimal basis. The current value of each objective coefficient is
     also displayed for reference.

     - *objective coefficient sensitivity* {flux : (reduced_cost, lower_obj_sensitivity, coeff_value, upper_obj_sensitivity)}
     - *rhs sensitivity* {constraint : (low, value, high)}
     - *bound sensitivity ranges* {flux : (lb_low, lb_high, ub_low, ub_high)}

    """
    SENSE_RHS = {}
    SENSE_BND = {}
    SENSE_OBJ = {}
    c_names = c.linear_constraints.get_names()
    rhs_val = c.linear_constraints.get_rhs()
    j_names = c.variables.get_names()

    rhs_sense = c.solution.sensitivity.rhs()
    bnd_sense = c.solution.sensitivity.bounds()
    obj_sense = c.solution.sensitivity.objective()
    obj_coeff = c.objective.get_linear()
    red_cost = c.solution.get_reduced_costs()

    for r in range(c.variables.get_num()):
        SENSE_BND.update({j_names[r] : (bnd_sense[r][0], bnd_sense[r][1], bnd_sense[r][2], bnd_sense[r][3])})
        SENSE_OBJ.update({j_names[r] : (red_cost[r], obj_sense[r][0], obj_coeff[r], obj_sense[r][1])})

    for s in range(c.linear_constraints.get_num()):
        SENSE_RHS.update({c_names[s] : (rhs_sense[s][0], rhs_val[s], rhs_sense[s][1])})

    return (SENSE_OBJ, SENSE_RHS, SENSE_BND)

def cplx_getDualValues(c):
    """
    Get the get the dual values of the solution

     - *c* a CPLEX LP

    Output is a dictionary of {name : value} pairs

    """
    d_names = c.linear_constraints.get_names()
    d_values = c.solution.get_dual_values()
    output = {}
    for j in range(len(d_names)):
        output.update({d_names[j] : d_values[j]})
    return output


def cplx_getSolutionStatus(c):
    """
    Returns one of:

     - *LPS_OPT*: solution is optimal;
     - *LPS_FEAS*: solution is feasible;
     - *LPS_INFEAS*: solution is infeasible;
     - *LPS_NOFEAS*: problem has no feasible solution;
     - *LPS_UNBND*: problem has unbounded solution;
     - *LPS_UNDEF*: solution is undefined.
     - *LPS_NONE*: no solution

    """
    # solution.get_status() returns an integer code
    if c.solution.get_solution_type() == c.solution.type.none:
        print("\nNo solution available\n")
        return 'LPS_NONE'
    status = c.solution.get_status()
    if status == c.solution.status.optimal:
        if not CPLX_SILENT_MODE:
            print("Model is optimal")
        return 'LPS_OPT'
    elif status == c.solution.status.feasible:
        print("\nModel is feasible")
        return 'LPS_FEAS'
    elif status == c.solution.status.unbounded:
        print("\nModel is unbounded")
        return 'LPS_UNBND'
    elif status == c.solution.status.infeasible:
        print("\nModel is infeasible")
        return 'LPS_INFEAS'
    elif status == c.solution.status.infeasible_or_unbounded:
        print("\nModel is infeasible or unbounded")
        return 'LPS_INFEAS or LPS_UNBND'
    elif status == c.solution.status.MIP_optimal:
        print('MILP optimal')
        return 'MILP_OPT'
    elif status == c.solution.status.MIP_optimal:
        print('MILP optimal')
        return 'MILP_OPT'
    elif status == c.solution.status.optimal_tolerance:
        print('MILP optimal within gap tolerance')
        return 'MILP_OPTTOL'
    else:
        return 'LPS_NONE'

def cplx_setObjective(c, pid, expr=None, sense='maximize', reset=True):
    """
    Set a new objective function note that there is a major memory leak in
    `c.variables.get_names()` whch is used when reset=True. If this is a problem
    use cplx_setObjective2 which takes *names* as an input:

     - *c* a CPLEX LP object
     - *pid* the r_id of the flux to be optimized
     - *expr* a list of (coefficient, flux) pairs
     - *sense* 'maximize'/'minimize'
     - *reset* [default=True] reset all objective function coefficients to zero

    """
    sense = sense.lower()
    if sense == 'max': sense = 'maximize'
    if sense == 'min': sense = 'minimize'
    if sense in ['maximise', 'minimise']:
        sense = sense.replace('se', 'ze')
    assert sense in ['maximize', 'minimize'], "\nsense must be ['maximize', 'minimize'] not %s" % sense
    c.objective.set_name(pid)
    if reset:
        vartemp = c.variables.get_names()
        new_obj = []
        for e in range(c.variables.get_num()):
            new_obj.append((vartemp[e], 0.0))
        c.objective.set_linear(new_obj)
        del new_obj, vartemp
    if expr != None:
        expr = [(e[1], e[0]) for e in expr]
        if __DEBUG__: print(expr)
        #print expr
        c.objective.set_linear(expr)
        del expr

    if sense == 'minimize':
        c.objective.set_sense(c.objective.sense.minimize)
        if __DEBUG__: print('Set minimizing')
    else:
        c.objective.set_sense(c.objective.sense.maximize)
        if __DEBUG__: print('Set maximizing')

def cplx_setObjective2(c, pid, names, expr=None, sense='maximize', reset=True):
    """
    Set a new objective function. This is a workaround function to avoid the
    e is a major memory leak in `c.variables.get_names()` whch is used
    in cplx_setObjective()  when reset=True.

    """
    sense = sense.lower()
    if sense == 'max': sense = 'maximize'
    if sense == 'min': sense = 'minimize'
    if sense in ['maximise', 'minimise']:
        sense = sense.replace('se', 'ze')
    assert sense in ['maximize', 'minimize'], "\nsense must be ['maximize', 'minimize'] not %s" % sense
    c.objective.set_name(pid)
    if reset:
        ##  var_num = c.variables.get_num()
        ##  names = range(var_num)
        new_obj = []
        for e in range(len(names)):
            new_obj.append((names[e], 0.0))
        c.objective.set_linear(new_obj)
        del new_obj, names
    if expr != None:
        expr = [(e[1], e[0]) for e in expr]
        if __DEBUG__: print(expr)
        c.objective.set_linear(expr)
        del expr
    if sense == 'minimize':
        c.objective.set_sense(c.objective.sense.minimize)
        if __DEBUG__: print('Set minimizing')
    else:
        c.objective.set_sense(c.objective.sense.maximize)
        if __DEBUG__: print('Set maximizing')


def cplx_setSingleConstraint(c, pid, expr, sense='E', rhs=0.0):
    ind = []
    val = []
    for e in expr:
        ind.append(e[1])
        val.append(e[0])
    lexp = cplex.SparsePair(ind, val)
    c.linear_constraints.add(lin_expr=[lexp], senses=sense, rhs=[rhs], names=[pid])
    del ind, val, lexp

def cplx_setBounds(c, pid, lb=None, ub=None):
    if lb != None:
        c.variables.set_lower_bounds(pid, lb)
    if ub != None:
        c.variables.set_upper_bounds(pid, ub)

def cplx_writeLPtoLPTfile(c, filename, title=None, Dir=None):
    """
    Write out a CPLEX model as an LP format file

    """
    if Dir != None:
        filename = os.path.join(Dir, filename)
    if title != None:
        c.set_problem_name(title)
    c.write(filename+'.lp', filetype='lp')
    print('LP output as {}'.format(filename+'.lp'))

def cplx_runInputScan(fba, exDict, wDir, input_lb=-10.0, input_ub=0.0, writeHformat=False, rationalLPout=False):
    """
    scans all inputs

    """
    debug_dir = os.path.join(wDir, 'lpt')
    ine_dir = os.path.join(wDir, 'ine')
    rat_dir = os.path.join(wDir, 'eslv')
    if not os.path.exists(debug_dir):
        os.mkdir(debug_dir)
    if writeHformat and not os.path.exists(ine_dir):
        os.mkdir(ine_dir)
    if rationalLPout and not os.path.exists(rat_dir):
        os.mkdir(rat_dir)

    optimal_growth_rates = {}
    infeasible_inputs = []
    modname = fba.sourcefile
    fbaid0 = fba.id
    fba.buildStoichMatrix()
    for inpt in exDict:
        fba.id = '%s(%s)' % (modname, inpt)
        ilb = fba.getFluxBoundByReactionID(inpt, 'lower').value
        iub = fba.getFluxBoundByReactionID(inpt, 'upper').value
        fba.setReactionBound(inpt, input_lb, 'lower')
        fba.setReactionBound(inpt, input_ub, 'upper')
        #fbalp = cplx_getModelFromObj(fba)
        fbalp = cplx_constructLPfromFBA(fba, fname=None)
        cplx_Solve(fbalp)
        if cplx_getSolutionStatus(fbalp) == 'LPS_OPT':
            opt = cplx_getOptimalSolution(fbalp)[2]
            if __DEBUG__: print('Optimal growth rate %s(%s) = %s' % (modname, inpt, opt))
            optimal_growth_rates.update({inpt : opt})
            if writeHformat:
                CBWrite.WriteModelHFormatFBA(fba, ine_dir, use_rational=False)
                CBWrite.WriteModelHFormatFBA(fba, ine_dir, use_rational=True)
            if rationalLPout:
                tmp = CBWrite.writeModelLP(fba, work_dir=rat_dir)
        else:
            print('CPLX returned an error code: {}'.format(cplx_getSolutionStatus(fbalp)))
            infeasible_inputs.append(inpt)
        fba.setReactionBound(inpt, ilb, 'lower')
        fba.setReactionBound(inpt, iub, 'upper')
        del fbalp
    fba.id = fbaid0
    return optimal_growth_rates, infeasible_inputs

##'''
##def cplx_runInputScanWithFVA(fba, exDict, wDir, input_lb=-10.0, input_ub=0.0, writeHformat=False, rationalLPout=False, NaNonFail=False):
    ##"""
    ##scans all inputs

    ##"""
    ##raise NotImplementedError('Not yet anyway')

    ##debug_dir = os.path.join(wDir, 'lpt')
    ##ine_dir = os.path.join(wDir, 'ine')
    ##rat_dir = os.path.join(wDir, 'eslv')
    ##if not os.path.exists(debug_dir):
        ##os.mkdir(debug_dir)
    ##if writeHformat and not os.path.exists(ine_dir):
        ##os.mkdir(ine_dir)
    ##if rationalLPout and not os.path.exists(rat_dir):
        ##os.mkdir(rat_dir)

    ##optimal_growth_rates = {}
    ##infeasible_inputs = []
    ##modname = fba.sourcefile
    ##fbaid0 = fba.id
    ##fva_out = []
    ##for inpt in exDict:
        ##fba.id = '%s(%s)' % (modname, inpt)
        ##ilb = fba.getFluxBoundByReactionID(inpt, 'lower').value
        ##iub = fba.getFluxBoundByReactionID(inpt, 'upper').value
        ##fba.setReactionBound(inpt, input_lb, 'lower')
        ##fba.setReactionBound(inpt, input_ub, 'upper')
        ##LPfile = CBWrite.writeModelLP(fba, work_dir=debug_dir)
        ##fbalp = PyscesGLPK.glpk_getModelFromLPT(LPfile, debug_dir)
        ####  fbalp.solve()
        ##PyscesGLPK.glpk_Solve(fbalp)
        ##if fbalp.solution_status == 'LPS_OPT':
            ##opt = PyscesGLPK.glpk_getOptimalSolution(fbalp)[2]
            ##if __DEBUG__: print 'Optimal growth rate %s(%s) = %s' % (modname, inpt, opt)
            ##optimal_growth_rates.update({inpt : opt})
            ##if writeHformat:
                ##CBWrite.WriteModelHFormatFBA(fba, ine_dir, realnum=False)
                ##CBWrite.WriteModelHFormatFBA(fba, ine_dir, realnum=True)
            ##if rationalLPout:
                ##tmp = CBWrite.writeModelLP(fba, work_dir=rat_dir)
            ##PyscesGLPK.glpk_addOptimalSolutionAsConstraint(fbalp, tol=None, bound='lower', optPercentage=100.0, reset_of=True)
            ##OUTPUT_ARRAY, ncol_names = doFluxVariability(fbalp, ignore_fluxes=[], write_debug_lps=False,\
                                        ##pause_on_skip=True, model_name=fba.id, work_dir=wDir)
            ##fva_out.append((inpt,OUTPUT_ARRAY, ncol_names))
            ##del OUTPUT_ARRAY, ncol_names
        ##else:
            ##print 'GLPK returned an error code: %s' % fbalp.solution_status
            ##infeasible_inputs.append(inpt)
            ##os.remove(LPfile)
        ##fba.setReactionBound(inpt, ilb, 'lower')
        ##fba.setReactionBound(inpt, iub, 'upper')
    ##fba.id = fbaid0
    ##return optimal_growth_rates, infeasible_inputs, fva_out
##'''

def cplx_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=False, oldlpgen=False, with_reduced_costs='unscaled', method='o'):
    """
    This is a utility function that does a presolve for FVA, MSAF etc. Generates properly formatted
    empty objects if pre_opt == False

     - *pre_opt* a boolean
     - *fba* a CBModel object
     - *objF2constr* add objective function as constraint
     - *quiet* [default=False] supress cplex output
     - *with_reduced_costs* [default='unscaled'] can be 'scaled' or 'unscaled'
     - *method* [default='o'] choose the CPLEX method to use for solution, default is automatic. See CPLEX reference manual for details

       - 'o': auto
       - 'p': primal
       - 'd': dual
       - 'b': barrier (no crossover)
       - 'h': barrier
       - 's': sifting
       - 'c': concurrent


    Returns: pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS

    """
    if oldlpgen:
        cpx = cplx_getModelFromObj(fba)
    else:
        cpx = cplx_constructLPfromFBA(fba, fname=None)
    if quiet:
        cplx_setOutputStreams(cpx, mode=None)
    OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    if pre_opt:
        status = cplx_Solve(cpx, method=method)
        if cplx_getSolutionStatus(cpx) == 'LPS_OPT':
            print('Valid Presolution')
            OPTIMAL_PRESOLUTION = True
            pre_sol, pre_oid, pre_oval = cplx_getOptimalSolution(cpx)
            fba.objectives[fba.activeObjIdx].solution, fba.objectives[fba.activeObjIdx].value = pre_sol, pre_oval
            if with_reduced_costs == 'scaled':
                REDUCED_COSTS = cplx_getReducedCosts(cpx, scaled=True)
            elif with_reduced_costs == 'unscaled':
                REDUCED_COSTS = cplx_getReducedCosts(cpx, scaled=False)
        else:
            print('Invalid Presolution')
            OPTIMAL_PRESOLUTION = False
            pre_sol = {}
            for r in fba.reactions:
                pre_sol.update({r.getId() : 0.0})
                r.reduced_cost = 0.0
            pre_oval = 0.0
            pre_oid = 'None'
            raise RuntimeError('\nPresolve failed to optimize this model and cannot continue!')
    else:
        pre_sol = {}
        for r in fba.reactions:
            pre_sol.update({r.getId() : 0.0})
        if objF2constr:
            pre_oval = fba.objectives[fba.activeObjIdx].value
            pre_oid = fba.objectives[fba.activeObjIdx].getId()
        else:
            pre_oval = 0.0
            pre_oid = 'None'
    for r in fba.reactions:
        r.value = pre_sol[r.id]
    return cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS

def cplx_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, oval, tol, optPercentage):
    """
    Take the objective function and "optimum" value and add it as a constraint
     - *cpx* a cplex object
     - *oval* the objective value
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)

    """

    # generate new constraint from old objective value (use non-zero coefficients)
    LCS = cpx.objective.get_linear()
    LCN = cpx.variables.get_names()
    new_constraint = []
    for c in range(len(LCS)):
        if LCS[c] != 0.0:
            new_constraint.append((LCS[c], LCN[c]))

    if rhs_sense == 'equal':
        pass
    elif (cpx.objective.get_sense() == cpx.objective.sense.maximize) and (rhs_sense == 'upper'):
        print('\nWarning: RHS sense error: \"upper\" does not match \"maximize\" changing to \"lower\"')
        rhs_sense = 'lower'
        time.sleep(1)
    elif (cpx.objective.get_sense() == cpx.objective.sense.minimize) and (rhs_sense == 'lower'):
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
            ub = numpy.ceil(oval/tol)*tol*optPercentage/100.0
        else:
            ub = oval*(optPercentage/100.0)
        cplx_setSingleConstraint(cpx, 'ObjCstr', expr=new_constraint, sense='L', rhs=ub)
        ##  cplx_setSingleConstraint(cpx, 'ObjCstr', expr=[(1, pre_oid)], sense='L', rhs=ub)
    elif rhs_sense == 'lower':
        if tol != None:
            lb = numpy.floor(oval/tol)*tol*optPercentage/100.0
        else:
            lb = oval*(optPercentage/100.0)
        cplx_setSingleConstraint(cpx, 'ObjCstr', expr=new_constraint, sense='G', rhs=lb)
        ##  cplx_setSingleConstraint(cpx, 'ObjCstr', expr=[(1, pre_oid)], sense='G', rhs=lb)
    else:
        raise RuntimeError("\nInvalid RHS sense: %s" % rhs_sense)


def cplx_MinimizeSumOfAbsFluxes(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=False, debug=False, objective_coefficients=None, return_lp_obj=False, oldlpgen=False, with_reduced_costs=None, method='o'):
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
     - *objective_coefficients* [default=None] a dictionary of (reaction_id : float) pairs that provide the are introduced as objective coefficients to the absolute flux value. Note that the default value of the coefficient (non-specified) is +1.
     - *return_lp_obj* [default=False] off by default when enabled it returns the CPLEX LP object
     - *with_reduced_costs* [default=None] if not None should be 'scaled' or 'unscaled'
     - *method* [default='o'] choose the CPLEX method to use for solution, default is automatic. See CPLEX reference manual for details

       - 'o': auto
       - 'p': primal
       - 'd': dual
       - 'b': barrier (no crossover)
       - 'h': barrier
       - 's': sifting
       - 'c': concurrent


    With outputs:

     - *fba* an update instance of a CBModel. Note that the FBA model objective function value is the original value set as a constraint

    """

    if objective_coefficients == None:
        objective_coefficients = {}
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
        debug_dir = os.path.join(work_dir, 'DEBUG')
        if not os.path.exists(debug_dir):
            os.mkdir(debug_dir)

    # generate a presolution
    print('CPLEX is using solver option: "{}"'.format(method))

    cpx = OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = cplx_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, with_reduced_costs=with_reduced_costs, method=method)
    # if required add the objective function as a constraint
    if objF2constr:
        cplx_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, tol, optPercentage)

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
    ##  lbs = [0.0]*len(abs_selected_reactions)
    ##  ubs = [numpy.inf]*len(abs_selected_reactions)
    ##  cpx.variables.add(lb=lbs, ub=ubs, names=abs_selected_reactions)
    cpx.variables.add(names=abs_selected_reactions)

    #cpx.variables.add(obj=[], lb=[], ub=[], types='', names=[], columns=[])
    #obj_func = [(1, r) for r in abs_selected_reactions]

    ##J - abs_J <= 0
    ##J + abs_J >= 0

    lin_exprs = []
    senses = []
    names = []
    rhss = []
    range_values = []
    obj_func = []

    for r in range(len(selected_reactions)):
        name1 = abs_selected_reactions[r]+'_c1'
        lin_expr1 = cplex.SparsePair(ind=[selected_reactions[r], abs_selected_reactions[r]], val=[1.0, -1.0])
        sense1 = 'L'
        name2 = abs_selected_reactions[r]+'_c2'
        lin_expr2 = cplex.SparsePair(ind=[selected_reactions[r], abs_selected_reactions[r]], val=[1.0, 1.0])
        sense2 = 'G'
        names.append(name1)
        lin_exprs.append(lin_expr1)
        senses.append(sense1)
        rhss.append(0.0)
        range_values.append(0.0)
        names.append(name2)
        lin_exprs.append(lin_expr2)
        senses.append(sense2)
        rhss.append(0.0)
        range_values.append(0.0)
        ##  print selected_reactions[r], abs_selected_reactions[r]
        if selected_reactions[r] in objective_coefficients:
            obj_func.append((objective_coefficients[selected_reactions[r]], abs_selected_reactions[r]))
        else:
            obj_func.append((1, abs_selected_reactions[r]))

    cplx_setObjective(cpx, 'MAFS', obj_func, 'min', reset=True)
    cpx.linear_constraints.add(lin_expr=lin_exprs, senses=senses, rhs=rhss, range_values=range_values, names=names)
    if debug:
        cplx_writeLPtoLPTfile(cpx, 'MSAF_base_(%s)' % time.time(), title=None, Dir=debug_dir)
        ##  cplx_writeLPtoLPTfile(cpx, 'MSAF_base_%s' % time.time() , title=None, Dir=debug_dir)

    cplx_Solve(cpx, method=method)
    cplx_setFBAsolutionToModel(fba, cpx, with_reduced_costs=with_reduced_costs)
    cplx_setSolutionStatusToModel(fba, cpx)

    minSum = cpx.solution.get_objective_value()
    fba.setAnnotation('min_flux_sum', minSum)
    fba.getActiveObjective().setValue(STORED_OPT)
    print('\nMinimizeSumOfAbsFluxes objective value: {}\n'.format(minSum))
    if quiet:
        cplx_setOutputStreams(cpx, mode='default')
    if not return_lp_obj:
        return minSum
    else:
        return cpx


def cplx_setSolutionStatusToModel(m, lp):
    """
    Sets the lp solutions status to the CBMPy model

    """
    m.SOLUTION_STATUS_INT = lp.solution.get_status()
    m.SOLUTION_STATUS = cplx_getSolutionStatus(lp)

def cplx_MinimizeNumActiveFluxes(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=False, debug=False, objective_coefficients=None, return_lp_obj=False, populate=None, oldlpgen=False):
    """
    Minimize the sum of active fluxes, updates the model with the values of the solution and returns the value
    of the MILP objective function (not the model objective function which remains unchanged). If population mode is activated
    output is as described below:

     Min: sum(Bi)
      Bi = 0 -> Ci Ji = 0

     Such that:
      NJi = 0
      Jbio = opt

     where:
       Binary Bi

    Arguments:

     - *fba* an FBA model object
     - *selected reactions* [default=None] means use all reactions otherwise use the reactions listed here
     - *pre_opt* [default=True] attempt to presolve the FBA and report its results in the ouput, if this is diabled and *objF2constr* is True then the vid/value of the current active objective is used
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
       Note this does not necessarily mean the upper or lower bound, although practically it will. If in doubt use *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value * (optPercentage/100.0)
     - *work_dir* [default=None] the MSAF working directory for temporary files default = cwd+fva
     - *debug* [default=False] if True write out all the intermediate MSAF LP's into work_dir
     - *quiet* [default=False] if enabled supress CPLEX output
     - *objF2constr* [default=True] add the model objective function as a constraint using rhs_sense etc. If
       this is True with pre_opt=False then the id/value of the active objective is used to form the constraint
     - *objective_coefficients* [default=None] a dictionary of (reaction_id : float) pairs that provide the are introduced as objective coefficients to the absolute flux value. Note that the default value of the coefficient (non-specified) is +1.
     - *return_lp_obj* [default=False] off by default when enabled it returns the CPLEX LP object
     - *populate* [default=None] enable search algorithm to find multiple (sub)optimal solutions. Set with a tuple of (RELGAP=0.0, POPULATE_LIMIT=20, TIME_LIMIT=300) suggested values only.
       - *RELGAP* [default=0.0] relative gap to optimal solution
       - *POPULATE_LIMIT* [default=20] terminate when so many solutions have been found
       - *TIME_LIMIT* [default=300] terminate search after so many seconds important with higher values of *POPULATION_LIMIT*
     - *with_reduced_costs* [default='uncsaled'] can be 'scaled', 'unscaled' or anything else which is None

    With outputs:

     - *mincnt* the objective function value OR
     - *mincnt, cpx* the objective function and cplex model OR
     - *populate_data, mincnt* a population data set OR
     - *populate_data, mincnt, cpx* both the cps object and population data set

    depending on selected flags.

    """
    if objective_coefficients == None:
        objective_coefficients = {}
    if fba.SCALED_REDUCED_COSTS:
        with_reduced_costs = 'scaled'
    else:
        with_reduced_costs = 'unscaled'
    if work_dir == None:
        work_dir = os.getcwd()
    else:
        assert os.path.exists(work_dir), '\nWhat did you think would happen now!'
    if debug:
        debug_dir = os.path.join(work_dir, 'DEBUG')
        if not os.path.exists(debug_dir):
            os.mkdir(debug_dir)

    base_reaction_names = fba.getReactionIds()

    # generate a presolution
    cpx = OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = cplx_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, with_reduced_costs=with_reduced_costs)
    # if required add the objective function as a constraint
    if objF2constr:
        cplx_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, tol, optPercentage)

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

    # this removes the objective function (redundant I think)
    #fba_obj_ids = fba.getActiveObjective().getFluxObjectiveReactions()
    #print fba_obj_ids
    #for R in fba_obj_ids:
        #if R in selected_reactions:
            #selected_reactions.pop(selected_reactions.index(R))
    #del R, fba_obj_ids
    NUM_FLX = len(selected_reactions)
    print('Total number of reactions: {}'.format(NUM_FLX))
    bin_selected_reactions = ['bin_%s' % r for r in selected_reactions]
    #print bin_selected_reactions
    # bin_J in [0,1]
    ##  lbs = [0.0]*len(bin_selected_reactions)
    ##  ubs = [numpy.inf]*len(bin_selected_reactions)
    ##  cpx.variables.add(lb=lbs, ub=ubs, names=bin_selected_reactions)
    #cpx.variables.add(names=bin_selected_reactions)
    cpx.variables.add(names=bin_selected_reactions, lb=numpy.zeros(len(bin_selected_reactions)),\
                      ub=numpy.ones(len(bin_selected_reactions)))
    ##cpx.variables.add(obj=[], lb=[], ub=[], types='', names=[], columns=[])

    for bv in bin_selected_reactions:
        cpx.variables.set_types(bv, cpx.variables.type.binary)

    obj_func = []

    for r in range(len(selected_reactions)):
        name1 = 'ind_'+selected_reactions[r]
        ##  print selected_reactions[r]
        lin_expr1 = cplex.SparsePair(ind=[selected_reactions[r]], val=[1.0])
        sense1 = 'E'
        rhs1 = 0.0
        indvar1 = cpx.variables.get_indices(bin_selected_reactions[r])
        complemented1 = 1
        cpx.indicator_constraints.add(lin_expr=lin_expr1, sense=sense1, rhs=rhs1, indvar=indvar1, complemented=complemented1, name=name1)
        #obj_func.append((1, bin_selected_reactions[r]))

        if selected_reactions[r] in objective_coefficients:
            obj_func.append((objective_coefficients[selected_reactions[r]], bin_selected_reactions[r]))
        else:
            obj_func.append((1, bin_selected_reactions[r]))

    #print obj_func
    #cpx.write('dump.lp')
    cplx_setObjective(cpx, 'MNAF', obj_func, 'min', reset=True)

    if debug:
        cplx_writeLPtoLPTfile(cpx, 'MNAF_base_(%s)' % time.time(), title=None, Dir=debug_dir)
        ##  cplx_writeLPtoLPTfile(cpx, 'MNAV_base_%s' % time.time() , title=None, Dir=debug_dir)

    #cpx.write('dump.lp')

    if populate == None:
        cplx_SolveMILP(cpx) #  cpx.solve()
        cplx_setFBAsolutionToModel(fba, cpx, with_reduced_costs=None)
        setReducedCosts(fba, REDUCED_COSTS)
        minCnt = cpx.solution.get_objective_value()
        fba.getActiveObjective().setValue(STORED_OPT)
        if quiet:
            cplx_setOutputStreams(cpx, mode='default')
        print('\nMinimizeNumActiveFluxes objective value: {}'.format(minCnt))
        if not return_lp_obj:
            return round(minCnt, 2)
        else:
            return round(minCnt, 2), cpx
    else:
        RELGAP = populate[0] #0.0 # relative gap to optimal solution
        POPULATE_LIMIT = populate[1] #20
        TIME_LIMIT = populate[2] #300 # seconds
        INTENSITY = cpx.parameters.mip.pool.intensity.values.very_aggressive
        DIVERSITY = cpx.parameters.mip.pool.replace.values.diversity
        #DIVERSITY = cpx.parameters.mip.pool.replace.values.firstin_firstout
        #DIVERSITY = cpx.parameters.mip.pool.replace.values.worst_objective
        ABSGAP = 0.0
        cpx.parameters.mip.pool.relgap.set(RELGAP) # Gunnar
        cpx.parameters.mip.pool.absgap.set(ABSGAP) # check this
        cpx.parameters.mip.pool.intensity.set(INTENSITY) # get "all" (sub)optimal solutions
        cpx.parameters.mip.limits.populate.set(POPULATE_LIMIT)
        cpx.parameters.mip.pool.replace.set(DIVERSITY)
        cpx.parameters.timelimit.set(TIME_LIMIT)
        #cpx.solve()
        cpx.populate_solution_pool()

        population = []
        var_names = cpx.variables.get_names()
        var_num = cpx.variables.get_num()
        pop_names = cpx.solution.pool.get_names()
        pop_num = cpx.solution.pool.get_num()
        print('CPLEX solution pool: {}'.format(pop_num))

        #pop_bin = []
        for p in range(pop_num):
            sol = cpx.solution.pool.get_values(p)
            binSum = 0
            for j_ in range(len(var_names)-1,-1,-1):
                if var_names[j_] not in base_reaction_names:
                    a = sol.pop(j_)
                    var_names.pop(j_)
                    binSum += round(a, 2)
            #pop_bin.append(binSum)
            population.append(sol)
        population.insert(0, base_reaction_names)
        for j_ in range(len(population[0])):
            fba.getReaction(population[0][j_]).setValue(population[1][j_])
        try:
            setReducedCosts(fba, REDUCED_COSTS)
        except Exception as ex:
            print(ex)
        minCnt = cpx.solution.get_objective_value()
        print('\nMinimizeNumActiveFluxes objective value: {}\n'.format(round(minCnt, 2)))
        fba.getActiveObjective().setValue(STORED_OPT)
        if quiet:
            cplx_setOutputStreams(cpx, mode='default')
        if not return_lp_obj:
            print('\nINFO: population scan now returns population and objective')
            return population, round(minCnt, 2)
        else:
            return population, round(minCnt, 2), cpx


def cplx_FluxVariabilityAnalysis(fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=True, debug=False, oldlpgen=False, markupmodel=True, default_on_fail=False, roundoff_span=10, method='o'):
    """
    Perform a flux variability analysis on an fba model:

     - *fba* an FBA model object
     - *selected reactions* [default=None] means use all reactions otherwise use the reactions listed here
     - *pre_opt* [default=True] attempt to presolve the FBA and report its results in the ouput, if this is disabled and *objF2constr* is True then the rid/value of the current active objective is used
     - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)
     - *work_dir* [default=None] the FVA working directory for temporary files default = cwd+fva
     - *debug* [default=False] if True write out all the intermediate FVA LP's into work_dir
     - *quiet* [default=False] if enabled, supress CPLEX output
     - *objF2constr* [default=True] add the model objective function as a constraint using rhs_sense etc. If
       this is True with pre_opt=False then the id/value of the active objective is used to form the constraint
     - *markupmodel* [default=True] add the values returned by the fva to the reaction.fva_min and reaction.fva_max
     - *default_on_fail* [default=False] if *pre_opt* is enabled replace a failed minimum/maximum with the solution value
     - *roundoff_span* [default=10] number of digits is round off (not individual min/max values)
     - *method* [default='o'] choose the CPLEX method to use for solution, default is automatic. See CPLEX reference manual for details

       - 'o': auto
       - 'p': primal
       - 'd': dual
       - 'b': barrier (no crossover)
       - 'h': barrier
       - 's': sifting
       - 'c': concurrent


    Returns an array with columns: Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus and a list containing the row names.

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
    cpx = OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = cplx_func_GetCPXandPresolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, method=method)
    # if required add the objective function as a constraint
    if objF2constr:
        cplx_func_SetObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, tol, optPercentage)
    if debug:
        cplx_writeLPtoLPTfile(cpx, 'FVA_base', title=None, Dir=debug_dir)

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
    for Ridx in range(NUM_FLX):
        R = selected_reactions[Ridx]
        OUTPUT_NAMES.append(R)
        MIN_STAT = MAX_STAT = 0
        # MIN
        ##  cplx_setObjective(cpx, R, expr=None, sense='min', reset=True)
        # TODO: bgoli: see whether this also works with 'minimize'
        cplx_setObjective(cpx, 'min%s' % R, [(1, R)], 'min', reset=True)
        ##  cplx_setBounds(c, id, min=None, max=None) # think about this
        MIN_STAT = cplx_Solve(cpx, method=method)
        if MIN_STAT == 1: # solved
            min_oval = cpx.solution.get_objective_value()
        elif MIN_STAT == 2: # unbound
            min_oval = -numpy.Inf
        elif MIN_STAT == 3:
            #min_oval = pre_sol[R] # try this as infeasible means no solution outside optimum
            min_oval = numpy.NaN
        else: # other failure
            min_oval = numpy.NaN
        if debug:
            cplx_writeLPtoLPTfile(cpx, '%smin' % R, title='min%s=%s' % (R,min_oval), Dir=debug_dir)

        # MAX
        ##  cplx_setObjective(cpx, R, expr=None, sense='max', reset=True)
        cplx_setObjective(cpx, 'max%s' % R, expr=None, sense='max', reset=False)
        ##  cplx_setBounds(c, id, min=None, max=None) # think about this
        MAX_STAT = cplx_Solve(cpx, method=method)
        if MAX_STAT == 1: # solved
            max_oval = cpx.solution.get_objective_value()
        elif MAX_STAT == 2: # unbound
            max_oval = numpy.Inf
        elif MAX_STAT == 3: # infeasible
            #max_oval = pre_sol[R] # try this as infeasible means no solution outside optimum
            max_oval = numpy.NaN
        else: # other fail
            max_oval = numpy.NaN
        if debug:
            cplx_writeLPtoLPTfile(cpx, '%smax' % R, title='max%s=%s' % (R,max_oval), Dir=debug_dir)

        # enables using the default value as a solution if the solver fails
        if pre_opt and default_on_fail:
            if MAX_STAT > 1 and not MIN_STAT > 1:
                max_oval = pre_sol[R]
            if MIN_STAT > 1 and not MAX_STAT > 1:
                min_oval = pre_sol[R]

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
            #if MAX_STAT == 1:
                #REAC.fva_max = max_oval
            #else:
                #REAC.fva_max = None
            #if MIN_STAT == 1:
                #REAC.fva_min = min_oval
            #else:
                #REAC.fva_min = None
            if R in REDUCED_COSTS:
                REAC.reduced_costs = REDUCED_COSTS[R]
        if not quiet and MAX_STAT > 1 or MIN_STAT > 1:
            print('Solver fail for reaction \"{}\" (MIN_STAT: {} MAX_STAT: {})'.format(R, MIN_STAT, MAX_STAT))
        cntr += 1
        if cntr == 200:
            tcnt += cntr
            print('FVA has processed {} of {} reactions'.format(tcnt, NUM_FLX))
            cntr = 0
    if quiet:
        cplx_setOutputStreams(cpx, mode='default')
    del cpx
    print('\nSinglecore FVA took: {} min (1 process)\n'.format((time.time()-s2time)/60.))
    print('Output array has columns:')
    print('Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus')
    return OUTPUT_ARRAY, OUTPUT_NAMES



def cplx_MultiFluxVariabilityAnalysis(lp, selected_reactions=None, tol=1e-10, rhs_sense='lower', optPercentage=100.0, work_dir=None, debug=False):
    """
    Perform a flux variability analysis on a multistate LP

     - *lp* a multistate LP
     - *selected reactions* [default=None] means use all reactions otherwise use the reactions listed here
     - *pre_opt* [default=True] attempt to presolve the FBA and report its results in the ouput
     - *tol*  [default=1e-10] do floor/ceiling the objective function constraint, otherwise floor/ceil to *tol*
     - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
     - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)
     - *work_dir* [default=None] the FVA working directory for temporary files default = cwd+fva
     - *debug* [default=False] if True write out all the intermediate FVA LP's into work_dir
     - *bypass* [default=False] bypass everything and only run the min/max on lp

    and returns an array with columns::

     Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus

    and a list containing the row names.
    """

    gc.enable()
    pre_opt = True
    debug_dir = work_dir
    if debug_dir == None:
        debug_dir = os.path.join(os.getcwd(), 'fva')
    else:
        assert os.path.exists(debug_dir), '\nWhat did you think would happen now!'
        debug_dir = os.path.join(debug_dir, 'fva')
    if debug and not os.path.exists(debug_dir):
        os.mkdir(debug_dir)

    TFname = CBTools.createTempFileName()
    OPTIMAL_PRESOLUTION = None
    pre_sol = pre_oid = pre_oval = None
    REDUCED_COSTS = {}
    if pre_opt:
        cplx_SolveMILP(lp, auto_mipgap=True)
        status = cplx_getSolutionStatus(lp)
        print(status)
        assert lp.solution.get_status() == lp.solution.status.MIP_optimal, '\n Fail not CPXMIP_OPTIMAL'
        print('Valid FVA presolution')
        OPTIMAL_PRESOLUTION = True
        pre_sol, pre_oid, pre_oval = cplx_getOptimalSolution(lp)
        REDUCED_COSTS = {}

    print(pre_oid, pre_oval)


    # generate new constraint from old objective value (use non-zero coefficients)
    LCS = lp.objective.get_linear()
    LCN = lp.variables.get_names()

    #new_constraint = []
    #for c in range(len(LCS)):
        #if LCS[c] != 0.0:
            #if 'xvar_' not in  LCN[c] and 'absL_' not in  LCN[c]:
                ##print LCS[c], LCN[c]
                #new_constraint.append((LCS[c], LCN[c]))

    #print len(LCS), len(LCN), len(new_constraint)
    #lp.write(os.path.join(work_dir, TFname+'_DEBUG.1.lp'),'lp')


    ## set objective constraint
    #if rhs_sense == 'equal':
        #cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='E', rhs=pre_oval)
        ###  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='E', rhs=pre_oval)
    #elif rhs_sense == 'upper':
        #if tol != None:
            #ub =  numpy.ceil(pre_oval/tol)*tol*optPercentage/100.0
        #else:
            #ub = pre_oval*(optPercentage/100.0)
        #cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='L', rhs=ub)
        ###  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='L', rhs=ub)
    #elif rhs_sense == 'lower':
        #if tol != None:
            #lb =  numpy.floor(pre_oval/tol)*tol*optPercentage/100.0
        #else:
            #lb = pre_oval*(optPercentage/100.0)
        #cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='G', rhs=lb)
        ###  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='G', rhs=lb)

    #lp.write(os.path.join(work_dir, TFname+'_DEBUG.2.lp'),'lp')

    NUM_FLX = lp.variables.get_num()
    VARIABLE_NAMES = lp.variables.get_names()
    print('Number of variables: {}'.format(NUM_FLX))
    if selected_reactions != None:
        for r in selected_reactions:
            assert r in VARIABLE_NAMES, "\n%s is not a valid reaction name" % r
        NUM_FLX = len(selected_reactions)
        print('Number of user selected variables: {}'.format(NUM_FLX))
        OUTPUT_ARRAY = numpy.zeros((NUM_FLX, 5))

    OUTPUT_NAMES = []
    iter_cnt = 0

    #raw_input('\n*****\nMultiFVA ready, press <enter> to continue:')
    Rdump = file(os.path.join(work_dir, TFname+'_resTemp.txt'),'w')
    NEW_RUN = True
    ADD_MILP_CONSTR = False
    for Ridx in range(NUM_FLX):
        if selected_reactions != None:
            R = selected_reactions[Ridx]
        else:
            R = VARIABLE_NAMES[Ridx]

        if ADD_MILP_CONSTR:
            if not NEW_RUN:
                lp.linear_constraints.delete('ObjCstr')
            NEW_RUN = False
            new_constraint = []
            for c in range(len(LCS)):
                if LCS[c] != 0.0:
                    if 'xvar_' not in  LCN[c] and 'absL_' not in  LCN[c] and LCN[c] == 'zvar1'+R:
                        print('Skipping objective', R)
                    elif 'xvar_' not in  LCN[c] and 'absL_' not in  LCN[c]:
                        #print LCS[c], LCN[c]
                        new_constraint.append((LCS[c], LCN[c]))

            print(len(LCS), len(LCN), len(new_constraint))

            lp.write(os.path.join(work_dir, TFname+'_DEBUG.1.lp'),'lp')
            # set objective constraint
            if rhs_sense == 'equal':
                cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='E', rhs=pre_oval)
                ##  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='E', rhs=pre_oval)
            elif rhs_sense == 'upper':
                if tol != None:
                    ub =  numpy.ceil(pre_oval/tol)*tol*optPercentage/100.0
                else:
                    ub = pre_oval*(optPercentage/100.0)
                cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='L', rhs=ub)
                ##  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='L', rhs=ub)
            elif rhs_sense == 'lower':
                if tol != None:
                    lb =  numpy.floor(pre_oval/tol)*tol*optPercentage/100.0
                else:
                    lb = pre_oval*(optPercentage/100.0)
                cplx_setSingleConstraint(lp, 'ObjCstr', expr=new_constraint, sense='G', rhs=lb)
                ##  cplx_setSingleConstraint(lp, 'ObjCstr', expr=[(1, pre_oid)], sense='G', rhs=lb)

            lp.write(os.path.join(work_dir, TFname+'_DEBUG.2.lp'),'lp')
            print(new_constraint)
            print(R)


        Rdump.write('%s\t%s\t%2.3f\t' % (Ridx, R, pre_sol[R]))
        iter_cnt += 1
        OUTPUT_NAMES.append(R)
        MIN_STAT = MAX_STAT = 0
        # MIN
        ##  cplx_setObjective(lp, R, expr=None, sense='min', reset=True)
        cplx_setObjective2(lp, 'min%s' % R, VARIABLE_NAMES, [(1, R)], 'min', reset=True)
        ##  cplx_setBounds(c, id, min=None, max=None) # think about this

        lp.write(os.path.join(work_dir, TFname+'_DEBUG.3a.lp'),'lp')
        cplx_SolveMILP(lp)

        if lp.solution.get_status() == lp.solution.status.MIP_optimal:
            #DBG# min_sol, min_oid, min_oval = cplx_getOptimalSolution(lp)
            min_sol, min_oid, min_oval = cplx_getOptimalSolution2(lp, VARIABLE_NAMES)
            ##  min_sol, min_oid, min_oval = 0.0,0.0,0.0
            MIN_STAT = 1
        else:
            ##  min_sol, min_oid, min_oval = 0.0
            min_sol = None
            min_oid = None
            min_oval = numpy.NaN
        if debug:
            cplx_writeLPtoLPTfile(lp, TFname+'%smin' % R, title='min%s=%s' % (R,min_oval), Dir=debug_dir)
        Rdump.write('%2.3f\t' % round(min_oval,3))
        # MAX
        ##  cplx_setObjective(lp, 'max%s' % R, expr=None, sense='max', reset=False)
        cplx_setObjective2(lp, 'max%s' % R, VARIABLE_NAMES, expr=None, sense='max', reset=False)
        ##  cplx_setBounds(c, id, min=None, max=None) # think about this

        lp.write(os.path.join(work_dir, TFname+'_DEBUG.3b.lp'),'lp')
        cplx_SolveMILP(lp)

        if lp.solution.get_status() == lp.solution.status.MIP_optimal:
            #DBG# max_sol, max_oid, max_oval = cplx_getOptimalSolution(lp)
            max_sol, max_oid, max_oval = cplx_getOptimalSolution2(lp, VARIABLE_NAMES)
            ##  max_sol, max_oid, max_oval = 0.0,0.0,0.0
            MAX_STAT = 1
        else:
            ##  max_sol, max_oid, max_oval = 0.0
            max_sol = None
            max_oid = None
            max_oval = numpy.NaN
        if debug:
            cplx_writeLPtoLPTfile(lp, TFname+'%smax' % R, title='max%s=%s' % (R,max_oval), Dir=debug_dir)
        Rdump.write('%2.3f\t' % round(max_oval,3))
        Rdump.write('%2.3f\n' % round(abs(max_oval - min_oval),3))

        OUTPUT_ARRAY[Ridx,0] = pre_sol[R]
        OUTPUT_ARRAY[Ridx,1] = min_oval
        OUTPUT_ARRAY[Ridx,2] = max_oval
        OUTPUT_ARRAY[Ridx,3] = MIN_STAT
        OUTPUT_ARRAY[Ridx,4] = MAX_STAT
        del max_sol, max_oid, max_oval, min_sol, min_oid, min_oval
        if MIN_STAT == 0 or MIN_STAT == 0:
            print('Solver fail for reaction')
            print(R, OUTPUT_ARRAY[Ridx])
        if iter_cnt == 100:
            iter_cnt = 0
            ##  mdtempFile = os.path.join(debug_dir,'_memDump.lp')
            ##  lp.write(mdtempFile, filetype='lp')
            ##  del lp
            gc.collect()
            ##  reload(cplex)
            print('Collecting garbage')
            ##  time.sleep(2)
            ##  lp = cplex.Cplex(mdtempFile)
    Rdump.close()
    #for r in OUTPUT_ARRAY:
        #print r
    print('Output array has columns:')
    print('Reaction, Variability Min, Variability Max, MinStatus, MaxStatus')
    del lp
    return OUTPUT_ARRAY, OUTPUT_NAMES


def cplx_singleGeneScan(fba, r_off_low=0.0, r_off_upp=0.0, optrnd=8, altout=False):
    """
    Perform a single gene deletion scan

     - *fba* a model object
     - *r_off_low* the lower bound of a deactivated reaction
     - *r_off_upp* the upper bound of a deactivated reaction
     - *optrnd* [default=8] round off the optimal value
     - *altout* [default=False] by default return a list of gene:opt pairs, alternatively (True) return an extended result set including gene groups, optima and effect map

    """
    # cplex optimization
    if fba.__single_gene_effect_map__ == None:
        fba.createSingleGeneEffectMap()
    #lpx = cbm.analyzeModel(fba, return_lp_obj=True)
    lpx = cplx_constructLPfromFBA(fba)
    lpx.solve()
    wtOpt = lpx.solution.get_objective_value()

    Jmap = fba.__single_gene_effect_map__.pop('keyJ')
    Emap = list(fba.__single_gene_effect_map__)

    cplx_setOutputStreams(lpx, mode=None)

    base_names = lpx.variables.get_names()
    base_lower = lpx.variables.get_lower_bounds()
    base_upper = lpx.variables.get_upper_bounds()

    base_lower_bounds = []
    base_upper_bounds = []

    t_cplex = 0.0

    for n_ in range(len(base_names)):
        base_lower_bounds.append((base_names[n_], base_lower[n_]))
        base_upper_bounds.append((base_names[n_], base_upper[n_]))

    results = []

    for pr_ in range(len(Emap)):
        ## cplex optimized method
        t_x1 = time.time()
        lpx.variables.set_lower_bounds(base_lower_bounds)
        lpx.variables.set_upper_bounds(base_upper_bounds)
        new_lower = []
        new_upper = []
        for r_ in range(len(Emap[pr_])):
            if not Emap[pr_][r_]:
                new_lower.append((Jmap[r_], r_off_low))
                new_upper.append((Jmap[r_], r_off_upp))
        if len(new_lower) > 0:
            lpx.variables.set_lower_bounds(new_lower)
        if len(new_upper) > 0:
            lpx.variables.set_upper_bounds(new_upper)
        lpx.solve()
        if lpx.solution.get_status() == lpx.solution.status.optimal:
            results.append({ 'opt' : round(lpx.solution.get_objective_value(), optrnd),
                             'deletions' : fba.__single_gene_effect_map__[Emap[pr_]],
                             #'activities' : Emap[pr_]
                           })
        else:
            results.append({ 'opt' : float('nan'),
                             'deletions' : fba.__single_gene_effect_map__[Emap[pr_]],
                             #'activities' : Emap[pr_]
                           })
        t_cplex += time.time() - t_x1

    fba.__single_gene_effect_map__['keyJ'] = Jmap
    del lpx, fba
    print('\nSingle gene deletion scan: {} mins'.format(t_cplex/60.0))

    singleRes = []
    for d in results:
        for g in d['deletions']:
            if g not in singleRes:
                if g == 'wt':
                    lbl = 'N/A'
                elif numpy.isnan(d['opt']):
                    lbl = 'no-solution'
                elif abs(round(wtOpt, optrnd) - d['opt']) <= 1.0e-10:
                    lbl = 'silent'
                elif abs(d['opt']) <= 1.0e-10:
                    lbl = 'lethal'
                else:
                    lbl = 'partial'
                singleRes.append((g, d['opt'], lbl))
            else:
                print('\nINFO: duplicate gene Id detected and skipped')

    if not altout:
        return singleRes
    else:
        return results


def cplx_WriteFVAtoCSV(pid, fva, names, Dir=None, fbaObj=None):
    """
    Takes the resuls of a FluxVariabilityAnalysis method and writes it to a nice
    csv file. Note this method has been refactored to `CBWrite.WriteFVAtoCSV()`.

     - *pid* filename_base for the CSV output
     - *fva* FluxVariabilityAnalysis() OUTPUT_ARRAY
     - *names* FluxVariabilityAnalysis() OUTPUT_NAMES
     - *Dir* [default=None] if set the output directory for the csv files
     - *fbaObj* [default=None] if supplied adds extra model information into the output tables

    """

    print('\nThis method is now a wrapper for the newer CBWrite.WriteFVAtoCSV method, provided for back compatbility.\n')

    CBWrite.WriteFVAtoCSV(pid, fva, names, Dir, fbaObj)
