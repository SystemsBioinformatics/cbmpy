'''
Created on Nov 12, 2014

@author: arne
'''

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, gc, random
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
from cplex.exceptions import CplexError
from pyscescbm import CBWrite, CBTools
from pyscescbm.CBConfig import __CBCONFIG__ as __CBCONFIG__

from solver.CBSolver import SolverFactory, Solver


__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

"""
Sets the noise level of the solver CPLX_RESULT_STREAM can be::

 - *None* silent i.e. no output
 - *file* set solver to silent and output logs to *CPLX_RESULT_STREAM_FILE* cplex_output.log
 - *iostream* set solver to silent and output logs to *CPLX_RESULT_STREAM_IO* csio
 - *'default'* or anything else noisy with full output closes STREAM_IO and STREAM_FILE (default)
"""
global CPLX_RESULT_STREAM, CPLX_RESULT_STREAM_FILE, CPLX_RESULT_STREAM_IO, CPLX_SILENT_MODE
CPLX_RESULT_STREAM = 'default'
CPLX_RESULT_STREAM_FILE = None
CPLX_RESULT_STREAM_IO = None
CPLX_SILENT_MODE = True

class CPLEXSolver(Solver):
    def __init__(self, lp=None):
        if lp==None:
            self.lp = cplex.Cplex()
        else:
            self.lp = lp
    
    def fixConSense(self, operator):
        """
        Fixes the sense of inequality operators, returns corrected sense symbol
    
        - *operator* the operator to check
    
        """
    
        if operator in ['<=','<','L']:
            oPr = 'L'
        elif operator in ['>=','>','G']:
            oPr = 'G'
        elif operator in ['=','E']:
            oPr = 'E'
        else:
            raise RuntimeError('INFO: invalid operator: %s' % operator)
        return oPr

    
    def solve(self, method='o'):
        """ solves the current problem """
        alg = self.lp.parameters.lpmethod.values

        if method == "o":
            self.lp.parameters.lpmethod.set(alg.auto)
        elif method == "p":
            self.lp.parameters.lpmethod.set(alg.primal)
        elif method == "d":
            self.lp.parameters.lpmethod.set(alg.dual)
        elif method == "b":
            self.lp.parameters.lpmethod.set(alg.barrier)
            self.lp.parameters.barrier.crossover.set(self.lp.parameters.barrier.crossover.values.none)
        elif method == "h":
            self.lp.parameters.lpmethod.set(alg.barrier)
        elif method == "s":
            self.lp.parameters.lpmethod.set(alg.sifting)
        elif method == "c":
            self.lp.parameters.lpmethod.set(alg.concurrent)
        else:
            print("Unrecognized option, using automatic")
            self.lp.parameters.lpmethod.set(alg.auto)
        try:
            self.lp.solve()
        except cplex.exceptions.CplexSolverError as ex:
            print("Exception raised during solve:\n\t\"{}\"".format(ex))
            return None
        # solution.get_status() returns an integer code
        status = self.lp.solution.get_status()
        if status == self.lp.solution.status.optimal:
            if not CPLX_SILENT_MODE:
                print("INFO: Model is optimal:", status)
        elif status == self.lp.solution.status.unbounded:
            if not CPLX_SILENT_MODE:
                print("INFO: Model is unbounded:", status)
            return status
        elif status == self.lp.solution.status.infeasible:
            if not CPLX_SILENT_MODE:
                print("INFO: Model is infeasible:", status)
            return status
        elif status == self.lp.solution.status.infeasible_or_unbounded:
            if not CPLX_SILENT_MODE:
                print("INFO: Model is infeasible or unbounded:", status)
            return status
    
        if not CPLX_SILENT_MODE:
            s_method = self.lp.solution.get_method()
            print("Solution status = " , status, ":", end=" ")
            # the following line prints the status as a string
            print(self.lp.solution.status[status])
            print("Solution method = ", s_method, ":", end=" ")
            print(self.lp.solution.method[s_method])
    
        s_type = self.lp.solution.get_solution_type()
        if s_type == self.lp.solution.type.none:
            print("INFO: No solution available")
            return None
        else:
            if not CPLX_SILENT_MODE:
                print("Objective value = " , self.lp.solution.get_objective_value())
            return status
    
    def setMIPGapTolerance(self,tol):
        """
        Sets the the relative MIP gap tolerance
        """
        self.lp.parameters.mip.tolerances.mipgap.set(tol)
        print('New MIP gap tolerance', self.lp.parameters.mip.tolerances.mipgap.get())
    
    def setSolverParameters(self,parset):
        pass
    
    def solveMILP(self, auto_mipgap=False):
        """
        Solve and MILP
    
         - *auto_mipgap* auto decrease mipgap until mipgap == absmipgap
    
        """
        try:
            self.lp.solve()
            if self.lp.solution.get_status() == self.lp.solution.status.optimal_tolerance:
                if auto_mipgap:
                    while self.lp.parameters.mip.tolerances.mipgap.get() >= self.lp.parameters.mip.tolerances.absmipgap.get():
                        self.lp.parameters.mip.tolerances.mipgap.set(self.lp.parameters.mip.tolerances.mipgap.get()/10.0)
                        self.lp.solve()
                        if self.lp.solution.get_status() == self.lp.solution.status.MIP_optimal:
                            break
                    print('\n\nMILP solution gap tolerance set to: %s\n\n' % self.lp.parameters.mip.tolerances.mipgap.get())
                    time.sleep(5)
                if self.lp.solution.get_status() == self.lp.solution.status.optimal_tolerance:
                    print('\n\nMILP solution found with gap tolerance consider adjusting gap tolerance with cplx_setGapTolerance\n\n')
                    #raise RuntimeWarning
        except cplex.exceptions.CplexSolverError as ex:
            print("\nException raised during solve\n")
            print(ex)
    
    def getObjectiveValue(self):
        """ returns current objective value (typically valid after solve) """
        return self.lp.solution.get_objective_value()
    
    def getObjectiveId(self):
        """ returns the name of the current objective function """
        return self.lp.objective.get_name()
    
    def getSolutionStatus(self):
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
        if self.lp.solution.get_solution_type() == self.lp.solution.type.none:
            print("\nNo solution available\n")
            return 'LPS_NONE'
        status = self.lp.solution.get_status()
        if status == self.lp.solution.status.optimal:
            if not CPLX_SILENT_MODE:
                print("Model is optimal")
            return 'LPS_OPT'
        elif status == self.lp.solution.status.feasible:
            print("\nModel is feasible")
            return 'LPS_FEAS'
        elif status == self.lp.solution.status.unbounded:
            print("\nModel is unbounded")
            return 'LPS_UNBND'
        elif status == self.lp.solution.status.infeasible:
            print("\nModel is infeasible")
            return 'LPS_INFEAS'
        elif status == self.lp.solution.status.infeasible_or_unbounded:
            print("\nModel is infeasible or unbounded")
            return 'LPS_INFEAS or LPS_UNBND'
        elif status == self.lp.solution.status.MIP_optimal:
            print('MILP optimal')
            return 'MILP_OPT'
        elif status == self.lp.solution.status.MIP_optimal:
            print('MILP optimal')
            return 'MILP_OPT'
        elif status == self.lp.solution.status.optimal_tolerance:
            print('MILP optimal within gap tolerance')
            return 'MILP_OPTTOL'
        else:
            return 'LPS_NONE'

    def isDualFeasible(self):
        """ checks if problem has been solved to dual feasibility """
        pass
    
    def addLinearConstraint(self, name, coef, rhs, sense):
        """ adds an additional linear constraint.
        
        Warning: Adding constraints manually might be much slower than creating
        the whole model from a metabolic network at once
        """
        ind = []
        val = []
        for colname, colval in coef:
            ind.append(colname)
            val.append(colval)
        lexp = cplex.SparsePair(ind, val)
        self.lp.linear_constraints.add(lin_expr=[lexp], senses=sense, rhs=[rhs], names=[name])
        del ind, val, lexp
    
    def addVariables(self, names, lb=None, ub=None, obj=None):
        """ adds additional variables.
        
        - *names* list of variable names
        - *lb* list of lower bounds
        - *ub* list of upper bounds
        - *obj* list of objective values
        
        if only one variable should be added, the lists can also be replaced by
        plain values
        
        Warning: Adding constraints manually might be much slower than creating
        the whole model from a metabolic network at once
        """
        pass
    
    def setLowerBounds(self, bounds):
        """ sets lower bounds of given variables
        
        bounds should be a dictionary mapping variable names to the new bounds
        """
        lbounds = []
        for colName, colBound in bounds.items():
            lbounds.append((colName, colBound))
        self.lp.variables.set_lower_bounds(lbounds)
    
    def setUpperBounds(self, bounds):
        """ sets upper bounds of given variables
        
        bounds should be a dictionary mapping variable names to the new bounds
        """
        ubounds = []
        for colName, colBound in bounds.items():
            ubounds.append((colName, colBound))  
        self.lp.variables.set_upper_bounds(ubounds)
    
    def setObjective(self, name='obj', coef=None, sense='maximize', reset=True):
        """ sets the objective.
        
        - *name* name of the objective
        - *coef* dictionary mapping variable names to the new coefficients
        - *sense* objective sense, can be one of 'min' or 'max'
        - *reset* if all other objective coefficients should be reset
        
        Note that there is a major memory leak in
        `c.variables.get_names()` whch is used when reset=True. 
        A workaround is implemented in the old code of cplx_setObjective2 which 
        takes *names* as an input.
        """
        sense = sense.lower()
        if sense == 'max': sense = 'maximize'
        if sense == 'min': sense = 'minimize'
        if sense in ['maximise', 'minimise']:
            sense = sense.replace('se','ze')
        assert sense in ['maximize', 'minimize'], "\nsense must be ['maximize', 'minimize'] not %s" % sense
        self.lp.objective.set_name(name)
        if coef != None:
            if reset:
                variables = self.lp.variables.get_names()
                new_obj = []
                for varname in self.lp.variables:
                    if coef.hasKey(varname):
                        new_obj.append((varname, coef[varname]))
                    else:
                        new_obj.append((varname, 0.0))
                self.lp.objective.set_linear(new_obj)
                del new_obj, variables
            else:
                new_obj = []
                for varname, varcoef in coef.items():
                    new_obj.append((varname, varcoef))
        else:
            variables = self.lp.variables.get_names()
            new_obj = []
            for varname in self.lp.variables:
                new_obj.append((varname, 0.0))
            self.lp.objective.set_linear(new_obj)
            del new_obj, variables
    
        if sense == 'minimize':
            self.lp.objective.set_sense(self.lp.objective.sense.minimize)
            if __DEBUG__: print('Set minimizing')
        else:
            self.lp.objective.set_sense(self.lp.objective.sense.maximize)
            if __DEBUG__: print('Set maximizing')
            pass
    
    
    def getObjectiveCoef(self, obj=None):
        """ fetches the objective coefficients of the given variables 
        
        - obj: dictionary that has variable names as keys, the coefficients will
            be written in the corresponding values.
            If obj is None, a new dictionary for all variables is created
        """
        
    
    def setRHS(self, rhs):
        """ sets right-hand sides of given constraints
        
        rhs should be a dictionary mapping constraint names to the new rhs
        """
        pass
        
    def setSense(self, sense):
        """ sets the sense of given constraints
        
        sense should be a dictionary mapping constraint names to the new sense
        """
        pass

    def deleteVariables(self, variables):
        """ deletes the variables with specified names """
        pass
    
    def deleteConstraints(self, cons):
        """ delete constraints with specified names """
        pass
    
    def getSolution(self):
        """
        extract the primal solution
        """
        s_val = []
        s_name = []
        fba_sol = {}
        try:
            s_val = self.lp.solution.get_values()
            s_name = self.lp.variables.get_names()
            for n in range(len(s_name)):
                fba_sol.update({s_name[n] : s_val[n]})
        except Exception as ex:
            print(ex)
            print('WARNING: No solution to get')
            s_val = []
            s_name = []
            fba_sol = {}
        del s_name,s_val
        return fba_sol

    def getReducedCosts(self, scaled=False):
        """
        Extract ReducedCosts from LP and return as a dictionary 'Rid' : reduced cost
    
         - *scaled* scale the reduced cost by the optimal flux value
        """
        s_name = self.lp.variables.get_names()
        r_costs = self.lp.solution.get_reduced_costs()
        objf_val = self.lp.solution.get_objective_value()
        output = {}
        s_val = None
        if scaled:
            s_val = self.lp.solution.get_values()
        for v in range(len(s_name)):
            if scaled:
                try:
                    r_val = r_costs[v]*s_val[v]/objf_val
                except:
                    r_val = 0.0
            else:
                r_val = r_costs[v]
            output.update({s_name[v] : r_val})
        del s_name, r_costs, s_val
        return output
    
    def getShadowPrices(self):
        """
        Returns a dictionary of shadow prices containing 'Rid' : (lb, rhs, ub)
        """
        c_names = self.lp.linear_constraints.get_names()
        rhs_sense = self.lp.solution.sensitivity.rhs()
        rhs_val = self.lp.linear_constraints.get_rhs()
        output = {}
        for s in range(self.lp.linear_constraints.get_num()):
            output.update({c_names[s] : (rhs_sense[s][0], rhs_val[s], rhs_sense[s][1])})
        return output
    
    def getSensitivities(self):
        """
        Get the sensitivities of each constraint on the objective function with input
    
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
        c_names = self.lp.linear_constraints.get_names()
        rhs_val = self.lp.linear_constraints.get_rhs()
        j_names = self.lp.variables.get_names()
    
        rhs_sense = self.lp.solution.sensitivity.rhs()
        bnd_sense = self.lp.solution.sensitivity.bounds()
        obj_sense = self.lp.solution.sensitivity.objective()
        obj_coeff = self.lp.objective.get_linear()
        red_cost = self.lp.solution.get_reduced_costs()
    
        for r in range(self.lp.variables.get_num()):
            SENSE_BND.update({j_names[r] : (bnd_sense[r][0], bnd_sense[r][1], bnd_sense[r][2], bnd_sense[r][3])})
            SENSE_OBJ.update({j_names[r] : (red_cost[r], obj_sense[r][0], obj_coeff[r], obj_sense[r][1])})
    
        for s in range(self.lp.linear_constraints.get_num()):
            SENSE_RHS.update({c_names[s] : (rhs_sense[s][0], rhs_val[s], rhs_sense[s][1])})
    
        return (SENSE_OBJ, SENSE_RHS, SENSE_BND)
    
    def getDualValues(self):
        """
        Get the get the dual values of the solution
    
        Output is a dictionary of {name : value} pairs
    
        """
        d_names = self.lp.linear_constraints.get_names()
        d_values = self.lp.solution.get_dual_values()
        output = {}
        for j in range(len(d_names)):
            output.update({d_names[j] : d_values[j]})
        return output
    
    def write(self, filename, title=None, Dir=None):
        """ write problem to file"""
        if Dir != None:
            filename = os.path.join(Dir, filename)
        if title != None:
            self.lp.set_problem_name(title)
        self.lp.write(filename+'.lp', filetype='lp')
        print('LP output as {}'.format(filename+'.lp'))
    
    
    def setOutputStreams(self, mode='default'):
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
            self.lp.set_log_stream(None)
            self.lp.set_results_stream(None)
            self.lp.set_warning_stream(None)
            CPLX_SILENT_MODE = True
        elif CPLX_RESULT_STREAM == 'file':
            if CPLX_RESULT_STREAM_FILE == None:
                CPLX_RESULT_STREAM_FILE = file('cplex_output.log', 'w')
            self.lp.set_log_stream(CPLX_RESULT_STREAM_FILE)
            self.lp.set_results_stream(CPLX_RESULT_STREAM_FILE)
            CPLX_SILENT_MODE = True
        elif CPLX_RESULT_STREAM == 'iostream':
            if CPLX_RESULT_STREAM_IO == None:
                CPLX_RESULT_STREAM_IO = csio.StringIO()
            self.lp.set_log_stream(CPLX_RESULT_STREAM_IO)
            self.lp.set_results_stream(CPLX_RESULT_STREAM_IO)
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
            self.lp.set_log_stream(os.sys.stdout)
            self.lp.set_results_stream(os.sys.stdout)

class CPLEXFactory(SolverFactory):
    
    def createEmpty(self):
        """ creates a new empty solver instance """
        solver = CPLEXSolver()
        return solver
    
    def create(self, fba = None, fname=None):
        """creates a new solver instance.
        
        If an fba instance (metabolic network) is given, the solver is setup
        to solve the fba instance. 
        Otherwise an empty problem is returned for manual setup.
        - *fname* optional filename if defined writes out the constructed lp
        """
        
        if fba == None:
            solver = CPLEXSolver()
            return solver
        else:
            _Stime = time.time()
            # defines
            csense = 'E' # constraint sense
        
            # define model and add variables
            solver = CPLEXSolver()
            lp = solver.lp
            lp.set_problem_name('%s' % (fba.getPid()))
            lp.variables.add(names=fba.N.col)
            # define objective
            osense = fba.getActiveObjective().operation.lower()
            if osense in ['minimize', 'minimise', 'min']:
                lp.objective.set_sense(lp.objective.sense.minimize)
            elif osense in ['maximize', 'maximise', 'max']:
                lp.objective.set_sense(lp.objective.sense.maximize)
            else:
                raise RuntimeError('\n%s - is not a valid objective operation' % osense)
            lp.objective.set_name(fba.getActiveObjective().getPid())
            lp.objective.set_linear([(fo.reaction, fo.coefficient) for fo in fba.getActiveObjective().fluxObjectives])
        
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
            rhs =[]
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
                lin_expr.append(cplex.SparsePair([fba.N.col[c] for c in nz], [Nmat[n,c] for c in nz]))
                rhs.append(RHSmat[n])
                senses.append(solver.fixConSense(fba.N.operators[n]))
                names.append(fba.N.row[n])
            #print senses
            lp.linear_constraints.add(lin_expr=lin_expr, senses=senses, rhs=rhs, names=names)
            #print 'New style lc:', time.time() - tnew
        
            # add user defined constraints
            if fba.CM != None:
                ## the numpy way
                #t2new = time.time()
                lin_expr = []
                rhs =[]
                names = []
                senses = []
                for n in range(CMmat.shape[0]):
                    if not GOSCI:
                        nz = CMmat[n].nonzero()[0]
                    else:
                        nz = CMmat[n].nonzero()[1]
        
                    lin_expr.append(cplex.SparsePair([fba.CM.col[c] for c in nz], [CMmat[n,c] for c in nz]))
                    rhs.append(CMrhs[n])
                    senses.append(solver.fixConSense(fba.CM.operators[n]))
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
        
    def getModelFromLP(self, lptFile, Dir=None):
        """
        Load a LPT (CPLEX format) file and return a CPLX LP model
    
         - *lptfile* an CPLEX LP format file
         - *Dir* an optional directory
    
        """
        if Dir != None:
            assert os.path.exists(Dir), '\nIncorrect path'
            lptFile = os.path.join(Dir, lptFile)
        lp = cplex.Cplex(lptFile)
        return CPLEXSolver(lp)

    def getCPLEXModelFromLP(self, lptFile, Dir=None):
        """
        Load a LPT (CPLEX format) file and return a CPLX LP model
    
         - *lptfile* an CPLEX LP format file
         - *Dir* an optional directory
    
        """
        return self.getModelFromLP(lptFile, Dir)
    
    
    def getModelFromObj(self, fba):
        """
        Return a CPLEX object from a FBA model object (via LP file)
        """
        _Stime = time.time()
        modId = fba.id
        ##  fba.id = '_cplxtmp_%s_.tmp' % random.randint(0,5)
        fba.id = '_cplxtmp_.tmp'
        LPF = CBWrite.writeModelLP(fba, work_dir=os.getcwd())
        fba.id = modId
        x = self.getModelFromLP(LPF)
        print('\ncplx_getModelFromObj time: {}\n'.format(time.time() - _Stime))
        return x

    def minimizeNumActiveFluxes(self, fba, selected_reactions=None, pre_opt=True, tol=None, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=False, debug=False, objective_coefficients={}, return_lp_obj=False, populate=None, oldlpgen=False):
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
         - *objective_coefficients* [default={}] a dictionary of (reaction_id : float) pairs that provide the are introduced as objective coefficients to the absolute flux value. Note that the default value of the coefficient (non-specified) is +1.
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
        if fba.SCALED_REDUCED_COSTS:
            with_reduced_costs = 'scaled'
        else:
            with_reduced_costs = 'unscaled'
        if work_dir == None:
            work_dir = os.getcwd()
        else:
            assert os.path.exists(work_dir), '\nWhat did you think would happen now!'
        if debug:
            debug_dir = os.path.join(work_dir,'DEBUG')
            if not os.path.exists(debug_dir):
                os.mkdir(debug_dir)
    
        base_reaction_names = fba.getReactionIds()
    
        # generate a presolution
        cpx = OPTIMAL_PRESOLUTION = None
        pre_sol = pre_oid = pre_oval = None
        REDUCED_COSTS = {}
        cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = self.presolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, with_reduced_costs=with_reduced_costs)
        # if required add the objective function as a constraint
        if objF2constr:
            cpx.setObjectiveFunctionAsConstraint(rhs_sense, pre_oval, optPercentage)
    
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
        cpx.lp.variables.add(names=bin_selected_reactions,lb=numpy.zeros(len(bin_selected_reactions)),\
                          ub=numpy.ones(len(bin_selected_reactions)))
        """
        cpx.variables.add(obj=[], lb=[], ub=[], types='', names=[], columns=[])
        """
    
        for bv in bin_selected_reactions:
            cpx.lp.variables.set_types(bv, cpx.lp.variables.type.binary)
    
        obj_func = {}
    
        for r in range(len(selected_reactions)):
            name1 = 'ind_'+selected_reactions[r]
            ##  print selected_reactions[r]
            lin_expr1 = cplex.SparsePair(ind = [selected_reactions[r]], val = [1.0])
            sense1 = 'E'
            rhs1 = 0.0
            indvar1 = cpx.lp.variables.get_indices(bin_selected_reactions[r])
            complemented1 = 1
            cpx.lp.indicator_constraints.add(lin_expr=lin_expr1, sense=sense1, rhs=rhs1, indvar=indvar1, complemented=complemented1, name=name1)
            #obj_func.append((1, bin_selected_reactions[r]))
    
            if selected_reactions[r] in objective_coefficients:
                obj_func[bin_selected_reactions[r]] = objective_coefficients[selected_reactions[r]]
                #obj_func.append((objective_coefficients[selected_reactions[r]], bin_selected_reactions[r]))
            else:
                obj_func[bin_selected_reactions[r]] = 1
                #obj_func.append((1, bin_selected_reactions[r]))
    
        #print obj_func
        #cpx.write('dump.lp')
        cpx.setObjective('MNAF', obj_func, 'min', reset=True)
    
        if debug:
            cpx.writeLPtoLPTfile('MNAF_base_(%s)' % time.time() , title=None, Dir=debug_dir)
            ##  cplx_writeLPtoLPTfile(cpx, 'MNAV_base_%s' % time.time() , title=None, Dir=debug_dir)
    
        #cpx.write('dump.lp')
    
        if populate == None:
            cpx.solveMILP() #  cpx.solve()
            cpx.setFBAsolutionToModel(fba, with_reduced_costs=None)
            cpx.setReducedCostsToModel(fba, REDUCED_COSTS)
            minCnt = cpx.getObjectiveValue()
            fba.getActiveObjective().setValue(STORED_OPT)
            if quiet:
                cpx.setOutputStreams(mode='default')
            print('\nMinimizeNumActiveFluxes objective value: {}'.format(minCnt))
            if not return_lp_obj:
                return round(minCnt, 2)
            else:
                return round(minCnt, 2), cpx
        else:
            RELGAP = populate[0] #0.0 # relative gap to optimal solution
            POPULATE_LIMIT = populate[1] #20
            TIME_LIMIT =  populate[2] #300 # seconds
            INTENSITY = cpx.lp.parameters.mip.pool.intensity.values.very_aggressive
            DIVERSITY = cpx.lp.parameters.mip.pool.replace.values.diversity
            #DIVERSITY = cpx.parameters.mip.pool.replace.values.firstin_firstout
            #DIVERSITY = cpx.parameters.mip.pool.replace.values.worst_objective
            ABSGAP = 0.0
            cpx.lp.parameters.mip.pool.relgap.set(RELGAP) # Gunnar
            cpx.lp.parameters.mip.pool.absgap.set(ABSGAP) # check this
            cpx.lp.parameters.mip.pool.intensity.set(INTENSITY) # get "all" (sub)optimal solutions
            cpx.lp.parameters.mip.limits.populate.set(POPULATE_LIMIT)
            cpx.lp.parameters.mip.pool.replace.set(DIVERSITY)
            cpx.lp.parameters.timelimit.set(TIME_LIMIT)
            #cpx.solve()
            cpx.lp.populate_solution_pool()
    
            population = []
            var_names = cpx.lp.variables.get_names()
            var_num = cpx.lp.variables.get_num()
            pop_names = cpx.lp.solution.pool.get_names()
            pop_num = cpx.lp.solution.pool.get_num()
            print('CPLEX solution pool: {}'.format(pop_num))
    
            #pop_bin = []
            for p in range(pop_num):
                sol = cpx.lp.solution.pool.get_values(p)
                binSum = 0
                for j_ in range(len(var_names)-1,-1,-1):
                    if var_names[j_] not in base_reaction_names:
                        a = sol.pop(j_)
                        b = var_names.pop(j_)
                        binSum += round(a, 2)
                #pop_bin.append(binSum)
                population.append(sol)
            population.insert(0, base_reaction_names)
            for j_ in range(len(population[0])):
                fba.getReaction(population[0][j_]).setValue(population[1][j_])
            try:
                cpx.setReducedCosts(fba, REDUCED_COSTS)
            except Exception as ex:
                print(ex)
            minCnt = cpx.lp.solution.get_objective_value()
            print('\nMinimizeNumActiveFluxes objective value: {}\n'.format(round(minCnt, 2)))
            fba.getActiveObjective().setValue(STORED_OPT)
            if quiet:
                cpx.setOutputStreams(mode='default')
            if not return_lp_obj:
                print('\nINFO: population scan now returns population and objective')
                return population, round(minCnt, 2)
            else:
                return population, round(minCnt, 2), cpx