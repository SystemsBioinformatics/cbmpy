'''
Created on Nov 12, 2014

@author: arne
'''

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os
import time

import numpy

from cbmpy import CBWrite

class Solver():
    """ abstract solver stub to ease integration of new solvers """
    
    def copy(self):
        """ creates an independent copy of this solver"""
        raise NotImplementedError('abstract method')
    
    def solve(self):
        """ solves the current problem """
        pass
    
    def getObjectiveValue(self):
        """ returns current objective value (typically valid after solve) """
        pass
    
    def getObjectiveId(self):
        """ returns the name of the current objective function """
        pass
    
    def getSolutionStatus(self):
        """
        Returns one of:
    
         - *LPS_OPT*: solution is optimal;
         - *LPS_FEAS*: solution is feasible;
         - *LPS_INFEAS*: solution is infeasible;
         - *LPS_NOFEAS*: problem has no feasible solution;
         - *LPS_UNBND*: problem has unbounded solution;
         - *LPS_UNDEF*: solution is undefined.
    
        """
        pass
    
    def isFeasible(self):
        """ checks if problem has been solved to feasibility """
        status = self.getSolutionStatus();
        return (status == 'LPS_OPT') or (status == 'LPS_UNBD')
    
    def isOptimal(self):
        """ checks if problem has been solved to optimality 

        """
        return self.getSolutionStatus() == 'LPS_OPT'
    
    def isDualFeasible(self):
        """ checks if problem has been solved to dual feasibility """
        pass
    
    def addLinearConstraint(self, name, coef, rhs, sense):
        """ adds an additional linear constraint.
        
        Warning: Adding constraints manually might be much slower than creating
        the whole model from a metabolic network at once
        """
        pass
    
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
    
    def setBounds(self, bounds):
        """ sets lower bounds of given variables
        
        bounds should be a dictionary mapping variable names to the new bounds
        the new bounds for each variable are given by the pair (lb, ub)
        """
        pass
    
    def getBounds(self, bounds=None):
        """ fetches the bounds of the given variables 
        
        - bounds: dictionary that has variable names as keys, the bounds will
            be written in the corresponding values.
            If bounds is None, a new dictionary for all variables is created
        """
        raise NotImplementedError('subclasses must implement this')
    
    def setLowerBounds(self, bounds):
        """ sets lower bounds of given variables
        
        bounds should be a dictionary mapping variable names to the new bounds
        """
        pass
    
    def setUpperBounds(self, bounds):
        """ sets upper bounds of given variables
        
        bounds should be a dictionary mapping variable names to the new bounds
        """
        pass
    
    def setObjective(self, name='obj', coef=None, sense='maximize', reset=True):
        """ sets the objective
        
        - *name* name of the objective
        - *coef* dictionary mapping variable names to the new coefficients
        - *sense* objective sense, can be one of 'min' or 'max'
        - *reset* if all other objective coefficients should be reset
        """
        pass
    
    def getObjectiveCoef(self, obj=None):
        """ fetches the objective coefficients of the given variables 
        
        - obj: dictionary that has variable names as keys, the coefficients will
            be written in the corresponding values.
            If obj is None, a new dictionary for all variables is created
        """
        pass
    
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
        
        This returns a dictionary assigning to each variable its value
        """
        pass

    def getReducedCosts(self, scaled=False):
        """
        Extract ReducedCosts from LP and return as a dictionary 'Rid' : reduced cost
    
         - *scaled* scale the reduced cost by the optimal flux value
        """
        pass

    def getShadowPrices(self):
        """
        Returns a dictionary of shadow prices containing 'Rid' : (lb, rhs, ub)
        """
        pass
    
    def getDualValues(self):
        """
        Get the get the dual values of the solution
    
        Output is a dictionary of {name : value} pairs
    
        """
        pass
    
    def setSolutionToModel(self, fba, with_reduced_costs='unscaled'):
        """
        Sets the FBA solution from a CPLEX solution to an FBA object
    
         - *fba* and fba object
         - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
           or anything else which is interpreted as None. Scaled is: s_rcost = (r.reduced_cost*rval)/obj_value
    
        """
        sol = self.getSolution()
        objval = self.getObjectiveValue()
        if self.isFeasible():
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
            RC = self.getReducedCosts(scaled=scaled)
            self.setReducedCostsToModel(fba, RC)
        else:
            self.setReducedCostsToModel(fba, {})

    def setReducedCostsToModel(self, fba, reduced_costs):
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
    
    def setObjectiveAsConstraint(self, rhs_sense, oval, optPercentage):
        """
        Take the objective function and "optimum" value and add it as a constraint
         - *cpx* a cplex object
         - *oval* the objective value
         - *tol*  [default=None] do not floor/ceiling the objective function constraint, otherwise round of to *tol* WEIRD PARAMETER, removed
         - *rhs_sense* [default='lower'] means objC >= objVal the inequality to use for the objective constraint can also be *upper* or *equal*
         - *optPercentage* [default=100.0] means the percentage optimal value to use for the RHS of the objective constraint: optimal_value*(optPercentage/100.0)
    
        """
    
        # generate new constraint from old objective value (use non-zero coefficients)
        con = self.getObjectiveCoef()
    
#         if rhs_sense == 'equal':
#             pass
#         elif cpx.obj.maximize and (rhs_sense == 'upper'):
#             print('\nWarning: RHS sense error: \"upper\" does not match \"maximize\" changing to \"lower\"')
#             rhs_sense = 'lower'
#             time.sleep(1)
#         elif not cpx.obj.maximize and (rhs_sense == 'lower'):
#             print('\nWarning: RHS sense error: \"lower\" does not match \"minimize\" changing to \"upper\"')
#             rhs_sense = 'upper'
#             time.sleep(1)
#         else:
#             print('\nRHS sense ok.')
    
        # set objective constraint
        if rhs_sense == 'equal':
            self.addLinearConstraint('ObjCstr', con, rhs=oval, sense='E')
        elif rhs_sense == 'upper':
#            if tol != None:
#                ub =  numpy.ceil(oval/tol)*tol*optPercentage/100.0
#            else:
            ub = oval*(optPercentage/100.0)
            self.addLinearConstraint('ObjCstr', con, rhs=ub, sense='L')
        elif rhs_sense == 'lower':
#            if tol != None:
#                lb =  numpy.floor(oval/tol)*tol*optPercentage/100.0
#            else:
            lb = oval*(optPercentage/100.0)
            self.addLinearConstraint(name='ObjCstr', coef=con, rhs=lb, sense='G')
        else:
            raise RuntimeError("\nInvalid RHS sense: %s" % rhs_sense)

    def write(self, filename, title=None, Dir=None):
        """ write problem to file"""
        pass
    
    def writeLPsolution(self, fba_sol, objf_name, objf_val, fname, Dir=None, separator=','):
        """
        This function writes the optimal solution, produced wth `cplx_getOptimalSolution` to file
    
         - *fba_sol* a dictionary of Flux : value pairs
         - *objf_name* the objective function flux id
         - *objf_val* the optimal value of the objective function
         - *fname* the output filename
         - *Dir* [default=None] use directory if not None
         - *separator* [default=','] the column separator
    
        """
        if Dir != None:
            assert os.path.exists(Dir), '\nPath does not exist'
            fname = os.path.join(Dir, fname)
        fname += '.csv'
        F = file(fname, 'w')
        cntr = 0
        F.write('%s%s%s\n' % ('ObjectiveFunction', separator, objf_name))
        for r in fba_sol:
            F.write('%s%s%f\n' % (r, separator, fba_sol[r]))
        F.flush()
        F.close()
        print('CSV exported to {}'.format(fname))

    def setOutputStreams(self, mode='default'):
        """
        Sets the noise level of the solver, mode can be one of:
    
         - *None* silent i.e. no output
         - *'file'* set solver to silent and output logs to some log file
         - *'iostream'* set solver to silent and output logs to some io stream
         - *'default'* or anything else noisy with full output closes STREAM_IO and STREAM_FILE (default)
        """
        pass

class SolverFactory():
    '''
    abstract solver factory stub with default implementations to ease 
    integration of new solvers.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        
    def createEmpty(self):
        """ creates a new empty solver instance """
        pass
    
    def create(self, fba = None, fname=None):
        """creates a new solver instance.
        
        If an fba instance (metabolic network) is given, the solver is setup
        to solve the fba instance. 
        Otherwise an empty problem is returned for manual setup.
        """
        pass

    def analyzeModel(self, fba, lpFname=None, return_lp_obj=False, with_reduced_costs='unscaled', with_sensitivity=False,\
                      build_n=True, quiet=False, method=None):
        """
        Optimize a model and add the result of the optimization to the model object
        (e.g. `reaction.value`, `objectiveFunction.value`). The stoichiometric
        matrix is automatically generated. This is a common function available in all
        solver interfaces. By default returns the objective function value
    
         - *fba* an instantiated PySCeSCBM model object
         - *lpFname* [default=None] the name of the intermediate LP file saved when this has a string value.
         - *return_lp_obj* [default=False] off by default when enabled it returns the PyGLPK LP object
         - *with_reduced_costs* [default='unscaled'] calculate and add reduced cost information to mode this can be: 'unscaled' or 'scaled'
           or anything else which is interpreted as 'None'. Scaled means s_rcost = (r.reduced_cost*rval)/obj_value
         - *with_sensitivity* [default=False] add solution sensitivity information (not yet implemented)
         - *del_intermediate* [default=False] delete the intermediary files after updating model object, useful for server applications
         - *build_n* [default=True] generate stoichiometry from the reaction network (reactions/reagents/species)
         - *quiet* [default=False] suppress solver output
         - *method* [default=solverDefault] select the solver method, see the respective solver documentation for details
    
        """
    
        if build_n:
            fba.buildStoichMatrix()
            #CBTools.addStoichToFBAModel(f)
        fid = fba.id
    
        if with_reduced_costs == 'scaled':
            fba.SCALED_REDUCED_COSTS = True
        elif with_reduced_costs == 'unscaled':
            fba.SCALED_REDUCED_COSTS = False
        else:
            fba.SCALED_REDUCED_COSTS = None
    
        if lpFname == None:
            flp = self.create(fba, fname=None)
            #f.id = '_glpktmp_.tmp'
        else:
            flp = self.create(fba, fname=lpFname)
            fid = lpFname
    
        _Stime = time.time()
    
        print('\nanalyzeModel FBA --> LP time: {}\n'.format(time.time() - _Stime))
        fba.id = fid
        if method != None:
            flp.solve(method)
        else:
            flp.solve()
    
        flp.setSolutionToModel(fba, with_reduced_costs=with_reduced_costs)
    
        fba.SOLUTION_STATUS =  flp.getSolutionStatus()
        # TODO: need to synchronise all solvers to same integer system
        if fba.SOLUTION_STATUS == 'LPS_OPT':
            fba.SOLUTION_STATUS_INT = 1
        else:
            fba.SOLUTION_STATUS_INT = 999
    
#        if oldlpgen and del_intermediate:
#            os.remove(LPF)
        objv = fba.getActiveObjective().getValue()
        print('\nanalyzeModel objective value: {}\n'.format(objv))
        if return_lp_obj:
            return flp
        else:
            del flp
            return objv
            
    def presolve(self, fba, pre_opt, objF2constr, quiet=False, oldlpgen=True, with_reduced_costs='unscaled'):
        """
        This is a utility function that does a presolve for FVA, MSAF etc. Generates properly formatted
        empty objects if pre_opt == False
    
         - *pre_opt* a boolean
         - *fba* a CBModel object
         - *objF2constr* add objective function as constraint
         - *quiet* [default=False] supress cplex output
         - *with_reduced_costs* [default='unscaled'] can be 'scaled' or 'unscaled'
    
        Returns: pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS
    
        """
    
        cpx = self.create(fba)
        OPTIMAL_PRESOLUTION = None
        pre_sol = pre_oid = pre_oval = None
        REDUCED_COSTS = {}
        if pre_opt:
            cpx.solve()
            if cpx.isOptimal():
                print('Valid Presolution')
                OPTIMAL_PRESOLUTION = True
                pre_oval = cpx.getObjectiveValue()
                pre_oid = cpx.getObjectiveId()
                pre_sol = cpx.getSolution()
                fba.objectives[fba.activeObjIdx].solution = pre_sol
                fba.objectives[fba.activeObjIdx].value = pre_oval
                if with_reduced_costs == 'scaled':
                    REDUCED_COSTS = cpx.getReducedCosts(scaled=True)
                elif with_reduced_costs == 'unscaled':
                    REDUCED_COSTS = cpx.getReducedCosts(scaled=False)
            else:
                print('Invalid Presolution')
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
    
    def fluxVariabilityAnalysis(self, fba, selected_reactions=None, pre_opt=True, objF2constr=True, rhs_sense='lower', optPercentage=100.0, work_dir=None, quiet=True, debug=False, oldlpgen=False, markupmodel=True, default_on_fail=False, roundoff_span=10):
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
        cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = self.presolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen)
        assert isinstance(cpx, Solver)
        # if required add the objective function as a constraint
        if objF2constr:
            cpx.setObjectiveFunctionAsConstraint(cpx, rhs_sense, pre_oval, optPercentage)
        if debug:
            cpx.writeLPtoLPTfile('FVA_base', title=None, Dir=debug_dir)
    
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
            cpx.setObjective('min%s' % R, {R:1}, 'min', reset=True)
            ##  cplx_setBounds(c, id, min=None, max=None) # think about this
            MIN_STAT = cpx.solve()
            if MIN_STAT == 1: # solved
                min_oval = cpx.getObjectiveValue()
            elif MIN_STAT == 2: # unbound
                min_oval = -numpy.Inf
            elif MIN_STAT == 3:
                #min_oval = pre_sol[R] # try this as infeasible means no solution outside optimum
                min_oval = numpy.NaN
            else: # other failure
                min_oval = numpy.NaN
            if debug:
                cpx.writeLPtoLPTfile('%smin' % R, title='min%s=%s' % (R,min_oval), Dir=debug_dir)
    
            # MAX
            ##  cplx_setObjective(cpx, R, expr=None, sense='max', reset=True)
            cpx.setObjective('max%s' % R, coef=None, sense='max', reset=False)
            ##  cplx_setBounds(c, id, min=None, max=None) # think about this
            MAX_STAT = cpx.solve()
            if MAX_STAT == 1: # solved
                max_oval = cpx.getObjectiveValue()
            elif MAX_STAT == 2: # unbound
                max_oval = numpy.Inf
            elif MAX_STAT == 3: # infeasible
                #max_oval = pre_sol[R] # try this as infeasible means no solution outside optimum
                max_oval = numpy.NaN
            else: # other fail
                max_oval = numpy.NaN
            if debug:
                cpx.writeLPtoLPTfile(cpx, '%smax' % R, title='max%s=%s' % (R,max_oval), Dir=debug_dir)
    
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
            cpx.setOutputStreams(mode='default')
        del cpx
        print('\nSinglecore FVA took: {} min (1 process)\n'.format((time.time()-s2time)/60.))
        print('Output array has columns:')
        print('Reaction, Reduced Costs, Variability Min, Variability Max, abs(Max-Min), MinStatus, MaxStatus')
        return OUTPUT_ARRAY, OUTPUT_NAMES
    
    def minimizeSumOfAbsFluxes(self, fba, selected_reactions=None, pre_opt=True, objF2constr=True, rhs_sense='lower', optPercentage=100.0, objective_coefficients={}, work_dir=None, quiet=False, debug=False, return_lp_obj=False, oldlpgen=False, with_reduced_costs=None):
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
        cpx = OPTIMAL_PRESOLUTION = None
        pre_sol = pre_oid = pre_oval = None
        REDUCED_COSTS = {}
        cpx, pre_sol, pre_oid, pre_oval, OPTIMAL_PRESOLUTION, REDUCED_COSTS = self.presolve(fba, pre_opt, objF2constr, quiet=quiet, oldlpgen=oldlpgen, with_reduced_costs=with_reduced_costs)
        assert isinstance(cpx, Solver)
        # if required add the objective function as a constraint
        if objF2constr:
            cpx.setObjectiveAsConstraint(cpx, rhs_sense, pre_oval, optPercentage)
    
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
        #basecols = len(cpx.cols)
        lb = NUM_FLX*[0]
        obj = NUM_FLX*[1]
        cpx.addVariables(names=abs_selected_reactions, lb=lb, obj=obj)
        # set objective name
        cpx.setObjective(name='MAFS', sense='min', reset=False)
        
        """
         J - abs_J <= 0
         J + abs_J >= 0
        """
        for r in range(NUM_FLX):
            conName = 'abs{}a'.format(selected_reactions[r])
            conCoef = {selected_reactions[r]:1,abs_selected_reactions[r]:-1}
            conRHS = 0
            conSense = '<='
            cpx.addLinearConstraint(conName, conCoef, conRHS, conSense)
            
            conName = 'abs{}b'.format(selected_reactions[r])
            conCoef = {selected_reactions[r]:1,abs_selected_reactions[r]:1}
            conRHS = 0
            conSense = '>='
            cpx.addLinearConstraint(conName, conCoef, conRHS, conSense)
                
    
        if debug:
            cpx.write(cpxlp=os.path.join(debug_dir, 'MSAF_base_(%s).lp' % time.time()))
            #cplx_writeLPtoLPTfile(cpx, , title=None, Dir=debug_dir)
            ##  cplx_writeLPtoLPTfile(cpx, 'MSAF_base_%s' % time.time() , title=None, Dir=debug_dir)
    
        cpx.solve()
        cpx.setSolutionToModel(fba, with_reduced_costs=with_reduced_costs)
    
        #cpx.write(cpxlp='test.lp'); return cpx
    
        minSum = cpx.getObjectiveValue()
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
        
    def runInputScan(self, fba, exDict, wDir, input_lb=-10.0, input_ub=0.0, writeHformat=False, rationalLPout=False):
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
        for input in exDict:
            fba.id = '%s(%s)' % (modname, input)
            ilb = fba.getFluxBoundByReactionID(input, 'lower').value
            iub = fba.getFluxBoundByReactionID(input, 'upper').value
            fba.setReactionBound(input, input_lb, 'lower')
            fba.setReactionBound(input, input_ub, 'upper')
            #fbalp = cplx_getModelFromObj(fba)
            fbalp = self.create(fba, fname=None)
            fbalp.solve()
            if fbalp.getSolutionStatus() == 'LPS_OPT':
                opt = fbalp.getObjectiveValue()
                if __DEBUG__: print('Optimal growth rate %s(%s) = %s' % (modname, input, opt))
                optimal_growth_rates.update({input : opt})
                if writeHformat:
                    CBWrite.WriteModelHFormatFBA(fba, ine_dir, use_rational=False)
                    CBWrite.WriteModelHFormatFBA(fba, ine_dir, use_rational=True)
                if rationalLPout:
                    tmp = CBWrite.writeModelLP(fba, work_dir=rat_dir)
            else:
                print('Solver returned an error code: {}'.format(fbalp.getSolutionStatus()))
                infeasible_inputs.append(input)
            fba.setReactionBound(input, ilb, 'lower')
            fba.setReactionBound(input, iub, 'upper')
            del fbalp
        fba.id = fbaid0
        return optimal_growth_rates, infeasible_inputs
        
# The following method is a bit more complicated to generalize
#     def singleGeneScan(self, fba, r_off_low=0.0, r_off_upp=0.0, optrnd=8, altout=False):
#         """
#         Perform a single gene deletion scan
#     
#          - *fba* a model object
#          - *r_off_low* the lower bound of a deactivated reaction
#          - *r_off_upp* the upper bound of a deactivated reaction
#          - *optrnd* [default=8] round off the optimal value
#          - *altout* [default=False] by default return a list of gene:opt pairs, alternatively (True) return an extended
#          result set including gene groups, optima and effect map
#     
#         """
#         # cplex optimization
#         if fba.__single_gene_effect_map__ == None:
#             fba.createSingleGeneEffectMap()
#         #lpx = cbm.analyzeModel(fba, return_lp_obj=True)
#         lpx = self.create(fba)
#         lpx.solve()
#         wtOpt = lpx.isOptimal()
#     
#         Jmap = fba.__single_gene_effect_map__.pop('keyJ')
#         Emap = list(fba.__single_gene_effect_map__)
#     
#         lpx.setOutputStreams(mode=None)
#     
#         base_names = lpx.variables.get_names()
#         base_lower = lpx.variables.get_lower_bounds()
#         base_upper = lpx.variables.get_upper_bounds()
#     
#         base_lower_bounds = []
#         base_upper_bounds = []
#     
#         t_cplex = 0.0
#     
#         for n_ in range(len(base_names)):
#             base_lower_bounds.append((base_names[n_], base_lower[n_]))
#             base_upper_bounds.append((base_names[n_], base_upper[n_]))
#     
#         results = []
#     
#         for pr_ in range(len(Emap)):
#             ## cplex optimized method
#             t_x1 = time.time()
#             lpx.variables.set_lower_bounds(base_lower_bounds)
#             lpx.variables.set_upper_bounds(base_upper_bounds)
#             new_lower = []
#             new_upper = []
#             for r_ in range(len(Emap[pr_])):
#                 if not Emap[pr_][r_]:
#                     new_lower.append((Jmap[r_], r_off_low))
#                     new_upper.append((Jmap[r_], r_off_upp))
#             if len(new_lower) > 0:
#                 lpx.variables.set_lower_bounds(new_lower)
#             if len(new_upper) > 0:
#                 lpx.variables.set_upper_bounds(new_upper)
#             lpx.solve()
#             if lpx.solution.get_status() == lpx.solution.status.optimal:
#                 results.append({ 'opt' : round(lpx.solution.get_objective_value(), optrnd),
#                                  'deletions' : fba.__single_gene_effect_map__[Emap[pr_]],
#                                  #'activities' : Emap[pr_]
#                                })
#             else:
#                 results.append({ 'opt' : float('nan'),
#                                  'deletions' : fba.__single_gene_effect_map__[Emap[pr_]],
#                                  #'activities' : Emap[pr_]
#                                })
#             t_cplex += time.time() - t_x1
#     
#         fba.__single_gene_effect_map__['keyJ'] = Jmap
#         del lpx, fba
#         print('\nSingle gene deletion scan: {} mins'.format(t_cplex/60.0))
#     
#         singleRes = []
#         for d in results:
#             for g in d['deletions']:
#                 if g not in singleRes:
#                     if numpy.isnan(d['opt']):
#                         lbl = 'no-solution'
#                     elif d['opt'] == wtOpt:
#                         lbl = 'silent'
#                     elif abs(d['opt']) <= 1.0e-11:
#                         lbl = 'lethal'
#                     else:
#                         lbl = 'partial'
#                     singleRes.append((g, d['opt'], lbl))
#                 else:
#                     print('\nINFO: duplicate gene Id detected and skipped')
#     
#         if not altout:
#             return singleRes
#         else:
#             return results
        
        
###############################################################################
# boilerplate for determining default solver
###############################################################################

from cbmpy.CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']
SOLVER_PREF = __CBCONFIG__['SOLVER_PREF']

__AVAILABLE_SOLVERS__ = []

__CBCONFIG__['SOLVER_ACTIVE'] = None
try:
    from . import CBGLPK
    __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK'
    __AVAILABLE_SOLVERS__.append('GLPK')
except ImportError:
    print('GLPK not available')
try:
    from . import CBCPLEX
    if __CBCONFIG__['SOLVER_ACTIVE'] == None:
        __CBCONFIG__['SOLVER_ACTIVE'] = 'CPLEX'
    else:
        __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK+CPLEX'
    __AVAILABLE_SOLVERS__.append('CPLEX')
except Exception as ex:
    print('\n\n')
    print(ex)
    print('\n\n')
    print('CPLEX not available')

if __CBCONFIG__['SOLVER_ACTIVE'] == None:
    print('\n*****\nWARNING: No linear solver present, please install IBM CPLEX with Python bindings or PyGLPK, please see http://cbmpy.sourceforge.net for Windows binary or http://tfinley.net/software/pyglpk for source.\n*****\n')

if __CBCONFIG__['SOLVER_ACTIVE'] == 'GLPK+CPLEX':
    if __CBCONFIG__['SOLVER_PREF'] == None:
        __CBCONFIG__['SOLVER_ACTIVE'] = None
        while __CBCONFIG__['SOLVER_ACTIVE'] not in ['CPLEX','GLPK']:
            __CBCONFIG__['SOLVER_ACTIVE'] = raw_input('\nSolver preference not set. Please select a solver [CPLEX or GLPK]: ')

    elif __CBCONFIG__['SOLVER_PREF'] == 'CPLEX':
        __CBCONFIG__['SOLVER_ACTIVE'] = 'CPLEX'
    elif __CBCONFIG__['SOLVER_PREF'] == 'GLPK':
        __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK'


###############################################################################
# old functions for compatibility 
###############################################################################

if __CBCONFIG__['SOLVER_ACTIVE'] == 'CPLEX':
    factory = CBCPLEX.CPLEXFactory()
elif __CBCONFIG__['SOLVER_ACTIVE'] == 'GLPK':
    factory = CBGLPK.GLPKFactory()
else:
    factory = None

def getFactory():
    """ fetches default factory """
    return factory

def createSolver(fba=None):
    """ creates a default solver. 
    
    *fba* metabolic network model to initialize the solver to solve FBA. If no
        FBA model is given, an empty LP isi initialized.
    """
    return factory.create(fba)

def analyzeModel(fba):
    factory.analyzeModel(fba)
    
def FluxVariabilityAnalysis(fba, selected_reactions=None, pre_opt=True, objF2constr=True, rhs_sense='lower', optPercentage=100.0):
    factory.fluxVariabilityAnalysis(fba, selected_reactions, pre_opt, objF2constr, rhs_sense, optPercentage)

def MinimizeSumOfAbsFluxes(fba, selected_reactions=None, pre_opt=True, objF2constr=True, rhs_sense='lower', optPercentage=100.0, objective_coefficients={}):
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

    With outputs:

     - *fba* an update instance of a CBModel. Note that the FBA model objective function value is the original value set as a constraint

    """
    factory.minimizeSumOfAbsFluxes(fba, selected_reactions, pre_opt, objF2constr, rhs_sense, optPercentage, objective_coefficients)