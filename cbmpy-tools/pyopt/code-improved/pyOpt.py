

import os, time, numpy, cplex, sys
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pyscescbm as cbm; import os



def readModel(file, maxDelete, bio_reaction, infinityValue, USE_GENE):
    """
    Function that reads in a model from SBML3 file
    Optional: function GrowthCondition.setCondition for setting flux values
        to mimic growth conditions
    
    - *file*                sbml filename
    - *bio_reaction*        inner objective, usually biomass reaction
    - *lp*                  primal LP, derived from model
    - *splitCouples*        list from reversible reactions that are split
    - *lpName*              name for LP, derived from filename
    - *GPRdict*             dictionary with (iso)enzymes mapped to reactions
    - *SubUdict*            dictionary with genes mapped to isoenzymes
    - *Scale*               upscale factor
    - *maxDelete*           maximum gene/reaction deletions allowed
    - *infinityValue*       boundary values that can be assumed infinite

    """

    # Read in model and set growth conditions
    model = cbm.CBRead.readSBML3FBC('Models/' + file)
    # remove equiliy bonds
    model.splitEqualityFluxBounds()
    model.setObjectiveFlux(bio_reaction, osense='maximize')
    lpName = file.replace('.xml', '')
    lpName = lpName.replace('SynthSinks/','')
    #model = GrowthCondition.setCondition(model, bio_reaction)
    cbm.analyzeModel(model)

    # Build dictionaries for GPR associations
    GPRdict = {}
    SubUdict = {}
    if USE_GENE:
        GPRdict, SubUdict= doGeneMapping(model)

    # Upscale reaction boundaries, if needed
    model, Scale = upscaleModel(model)
    infinityValue = infinityValue * Scale

    # Write LP
    cbm.CBWrite.writeModelLP(model, work_dir = 'LPs', fname = lpName + '_primal')
    lp = cplex.Cplex('LPs/' + lpName + '_primal.lp')

    return lp, lpName, GPRdict, SubUdict, Scale, maxDelete, infinityValue


def upscaleModel(model):
    """
    Function that scales a model in case of small reaction bounds
    
    This improves the accuracy of the solver, since small fluxes cannot
    be confused anymore with no flux
    
    - *model*               Pysces model
    - *Scale*               scale factor

    """

    # Get and map flux bounds for each reaction
    # Sort absolute values
    reactionBoundaries = []
    boundDict = {}
    for r_ in model.getReactionIds():
        lb = model.getReactionBounds(r_)[1]
        ub = model.getReactionBounds(r_)[2]
        reactionBoundaries.append(abs(lb))
        reactionBoundaries.append(abs(ub))
        boundDict[r_] = [lb, ub]
    reactionBoundaries = list(set(reactionBoundaries)- set([0.0]))
    reactionBoundaries = sorted(reactionBoundaries)

    # Determine the hight of the first value from the sorted values
    # if this value is too small, the model is scales by a factor of 10 or 100
    Scale = 1
    if reactionBoundaries[0] < 1:
        for r_ in boundDict:
            model.setReactionBounds(r_, boundDict[r_][0]*100, boundDict[r_][1]*100)
        Scale = 100

    elif reactionBoundaries[0] < 10:
        for r_ in boundDict:
            model.setReactionBounds(r_, boundDict[r_][0]*10, boundDict[r_][1]*10)
        Scale = 10

    return model, Scale


def doGeneMapping(model):
    """
    Function that maps enzymes and genes to reactions
    
    This function works only if the GPR associations are defined
    as follows: (g1 and g2 and g6) or ((g3 or g10) and g12)
    
    - *model*               Pysces model
    - *GPRdict*             dictionary with (iso)enzymes mapped to reactions
    - *SubUdict*            dictionary with genes mapped to isoenzymes
    
    """
    import re
    def unique_list(seq):
        """Function to remove duplicates"""
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    # Get GPR associations for reactions (strings) and
    # and split according to keywords 'and', 'or'
    GPRdict = {}
    model.createGeneAssociationsFromAnnotations()
    reactions = model.getReactionIds()
    SubUdict = {}
    no_associations = 0
    for r_ in reactions:
            try:
                ass = model.getReaction(r_).getAnnotations()
                if 'GENE ASSOCIATION' in ass:
                    g_ = ass['GENE ASSOCIATION']
                elif 'GENE_ASSOCIATION' in ass:
                    g_ = ass['GENE_ASSOCIATION']
                elif 'gene_association' in ass:
                    g_ = ass['gene_association']
                elif 'gene association' in ass:
                    g_ = ass['gene association']
                if g_ != 'None' and g_ != '' :
                    # Enzymes
                    # g_ = g_.split(') or (')
                    g_ = re.split(r'\)\s+or\s+\(|\)\s+or\s+|\s+or\s+\(',g_)
                    S_list = []
                    for enzyme in g_:
                        enzyme = enzyme.replace(')','')
                        enzyme = enzyme.replace('(','')
                        # Isoenzymes
                        enzyme = enzyme.replace(' or ','_or_')
                        # Subunits
                        subunits = enzyme.split(' and ')
                        subunits_mod = []
                        for s in subunits:
                            # remove extra space
                            tmp = s.replace(' ','')
                            # replace possible dashes
                            tmp = tmp.replace('-','_')
                            # add gene prefix
                            tmp = 'GENE_' + tmp
                            subunits_mod.append(tmp)

                        S_list.append(subunits_mod)
                
                    # Dictionary for isoenzymes
                    for enzymes in S_list:
                        for gene in enzymes:
                            gene = gene.replace(' ','')
                            if 'or' in gene:
                                # SubUdict[gene] = gene.split('_or_')
                                SubUdict[gene] = unique_list(gene.split('_or_'))
        
                # GPRdict[r_] = S_list
                GPRdict[r_] = [unique_list(s) for s in S_list]
            except:
                    no_associations+=1
    print '{} of {} reactions have no GPR Associations' .format(no_associations,len(reactions))
    # print GPRdict
    # print SubUdict
    # raw_input()

    return GPRdict, SubUdict


def solveInner(lp, objMinFactor):
    """
    Function that solves the inner objective of the bilevel lp
    Inner objective is usually maximization of biomass
    
    - *lp*                  primal lp object
    - *objMinFactor*        percentage of maximum inner objective that
                            should be respected by outer objective
    - *objVal*              'biomass' value, related to objMinFactor
    - *linCoeffVars*        variables of inner objective, usually one
                            (biomass reaction)

    """

    # calculate ObjFunc value and minimum required value
    lp.solve()
    objVal = lp.solution.get_objective_value()
    objVal = objMinFactor*objVal

    # define objective flux(es)
    linCoeffVars = lp.objective.get_linear()
    linCoeffVars = lp.variables.get_names([a for a in range(len(linCoeffVars)) if linCoeffVars[a] != 0.0])

    return objVal, linCoeffVars


def losenInnerObj(lp, infinityValue, linCoeffVars, objVal):
    """
    Modifies the flux bounds, such that the bilevel program becomes more 
    easily solvable, because flux boundaries lead to dual variables. 
    This approach losens the equality constraint of the primal and dual, 
    and therefore losens the inner optimization for biomass

    Sets the lower bound for the inner objective
    
    - *lp*                  primal lp object
    - *infinityValue*       bound values that can be considered infinite
    - *linCoeffVars*        inner objective variables
    - *objVal*              minumum value for inner objective
    - *fix_lower_bounds*    lower bounds that are added to the problem 
                            after the bilevel lp is formed
    - *fix_upper_bounds*    lower bounds that are added to the problem 
                            after the bilevel lp is formed
    
    """

    new_lower_bounds = []
    new_upper_bounds = []
    fix_lower_bounds = []
    fix_upper_bounds = []

    lplb = lp.variables.get_lower_bounds()
    lpub = lp.variables.get_upper_bounds()
    lpnam = lp.variables.get_names()

    # Make large bounds cplex.infinite
    for v_ in range(len(lpnam)):
        if lplb[v_] <= -infinityValue:
            new_lower_bounds.append((lpnam[v_], -cplex.infinity))

        if lpub[v_] >= infinityValue:
            new_upper_bounds.append((lpnam[v_], cplex.infinity))

        # Reversible reactions
        # if lplb[v_] < 0.0 and lpub[v_] > 0.0 :
        #    new_lower_bounds.append((lpnam[v_], -cplex.infinity))
        #    new_upper_bounds.append((lpnam[v_], cplex.infinity))
        #    if lplb[v_] > -infinityValue:
        #        fix_lower_bounds.append((lpnam[v_], lplb[v_]))

    # set new upper and lower bounds
    if new_lower_bounds:
        lp.variables.set_lower_bounds(new_lower_bounds)
    if new_upper_bounds:
        lp.variables.set_upper_bounds(new_upper_bounds)


    # set lower bound for inner objective in indicator method
    # this assumes that the inner objective is a maximization problem
    for v_ in linCoeffVars:
        lp.variables.set_lower_bounds(v_, objVal)

    del new_lower_bounds
    del lplb, lpnam, v_, lpub

    return lp, fix_lower_bounds, fix_upper_bounds


def characterizePrimal(lp):
    """
    Defines primal variables, and their boundaries within the 
    primal problem formulation (using function call). 
    This might enhance the speed, but is no requirement.
    
    - *lp*                  cplex lp object
    - *primalConstraints*   primal constraint names (metabolites)
    - *primalVariables*     all reactions

    """

    # *primalVariables* reactions
    # *primalConstraints* metabolite flux balances
    primalVariables = lp.variables.get_names()
    primalConstraints = lp.linear_constraints.get_names()


    # rename the primal constraints as p<name>
    for i_ in range(lp.linear_constraints.get_num()):
        lp.linear_constraints.set_names(i_, 'p'+lp.linear_constraints.get_names(i_))


    return primalConstraints, primalVariables


def cplx_buildDual(lp):
    """
    Build the dual problem and return as a cplex object

     - *lp*                 cplex lp object
     - *lpd*                dual of lp object

    """
    fname = '{}.mip'.format(str(time.time()).split('.')[0])
    lp.write(fname, 'dua')
    lpd = cplex.Cplex()
    lpd.read(fname)
    os.remove(os.path.abspath(fname))

    return lpd


def characterizeDual(lpd,primalConstraints, lpName):
    """
    Defines dual variables, and their boundaries within the
    primal problem formulation (using function call). 
    This might enhance the speed, but is no requirement.
    
    - *lpd*                 cplex dual lp object
    - *dualConstraints*     dual constraint names (reactions)
    - *dualVariables*       dual variables

    """

    dualVariables = []
    for v_ in lpd.variables.get_names():
        if v_ not in primalConstraints:
            dualVariables.append(v_)
    dualConstraints = lpd.linear_constraints.get_names()

    # rename the dual constraints as d<name>
    for i_ in range(lpd.linear_constraints.get_num()):
        lpd.linear_constraints.set_names(i_, 'd'+lpd.linear_constraints.get_names(i_))


    #write dual to file
    #lpd.write('03.lp')
    lpd.write(os.path.join('LPs', lpName+'_dual.lp'))

    return dualConstraints, dualVariables


def cplx_mergePrimalDual(lp, lpd, newobj=None, osense='maximize'):
    """
    Merges the primal and dual into a single new problem constrained with
    pdConst: obj(primal) = obj(dual)

     - *lp*                 primal lp object
     - *lpd*                dual lp object
     - *newobj*             the new combined objective function as a list of [('id', coefficient)] pairs
     - *osense*             the objective sense

    """
    # get primal and dual objective
    linCoeffVars = [x for x in zip(lp.variables.get_names(), lp.objective.get_linear()) if x[1] != 0.0]
    lpdobj = [x for x in zip(lpd.variables.get_names(), lpd.objective.get_linear()) if x[1] != 0.0]

    # build bounds and constraints from lp and lpd
    lpx = cplex.Cplex()
    lpx.variables.add(lb=lp.variables.get_lower_bounds(),\
                        ub=lp.variables.get_upper_bounds(), names=lp.variables.get_names())
    lpx.variables.add(lb=lpd.variables.get_lower_bounds(),\
                        ub=lpd.variables.get_upper_bounds(), names=lpd.variables.get_names())
    lpx.linear_constraints.add(lin_expr=lp.linear_constraints.get_rows(), rhs=lp.linear_constraints.get_rhs(),\
                               names=lp.linear_constraints.get_names(), senses=lp.linear_constraints.get_senses())
    newC = []
    for c_ in lpd.linear_constraints.get_rows():
        c_.ind = [lpx.variables.get_indices(x) for x in [lpd.variables.get_names(y) for y in c_.ind]]
        newC.append(c_)

    lpx.linear_constraints.add(lin_expr=newC, rhs=lpd.linear_constraints.get_rhs(),\
                               names=lpd.linear_constraints.get_names(), senses=lpd.linear_constraints.get_senses())

    # all variables
    cVar = lpx.variables.get_names()

    # Set constraint linking the primal and the dual objectives
    objCnstr = []
    for c in linCoeffVars:
        objCnstr.append((cVar.index(c[0]), c[1]))

    for c in lpdobj:
        objCnstr.append((cVar.index(c[0]), -c[1]))

    objCnstr = cplex.SparsePair(ind=[c[0] for c in objCnstr], val=[v[1] for v in objCnstr])
    lpx.linear_constraints.add(lin_expr=[objCnstr], rhs=[0.0], senses=['E'], names=['objLinkC'])


    # define objective of new bilevel MILP
    try:
        if newobj != None:
            newOname = 'pd_'
            for of_ in newobj:
                newOname += '{}_'.format(of_[0])
                lpx.objective.set_linear(of_[0], of_[1])
            lpx.objective.set_name(newOname[:-1])
        if osense == 'maximize':
            lpx.objective.set_sense(lpx.objective.sense.maximize)
        elif osense == 'minimize':
            lpx.objective.set_sense(lpx.objective.sense.minimize)
    except:
        sys.exit('\n\nThe bilevelObjective is not known. Make sure that the reaction name is correct and that the reaction is defined as <reactionID>_fwd or <reactionID>_rev, if using COBRA.\n\n')

    return lpx


def createBinariesReac(lppd, primalVariables, maxDelete):
    """
    Creates binaries for each primal variable (reaction)
    Sets a maximum constraint on reactions that can be knocked out
    
    - *lppd*            bilevel lp
    - *primalVariables* reactions
    - *maxDelete*       maximum deletions allowed
    """
    noBinReac = []
    lplb = lppd.variables.get_lower_bounds()
    lpub = lppd.variables.get_upper_bounds()
    lpnam = lppd.variables.get_names()

    # Negative upperbounds and positive lower bounds cannot be in fused primal objective,
    # therefore, they are not included to be knocked out.
    for v_ in range(len(lpnam)):
        if lplb[v_] > 0:
            noBinReac.append(lpnam[v_])

        if lpub[v_] < 0:
            noBinReac.append(lpnam[v_])

    binPrimal = []
    for v_ in primalVariables:
        if v_ not in noBinReac:
            lppd.variables.add([0], [0], [1], types='B', names=['bin_'+v_])
            binPrimal.append('bin_'+v_)

    # add number of minimal active fluxes or genes
    lppd.linear_constraints.add([cplex.SparsePair(binPrimal, numpy.ones(len(binPrimal)))],
                                 "G", rhs=[len(binPrimal)-maxDelete], names=['MinActR'])

    return binPrimal, lppd


def createBinariesGenes(lppd, GPRdict, SubUdict, maxDelete):
    """
    Creates binaries for each primal variable (reaction), isoenzyme and gene
    Sets a maximum constraint on genes that can be knocked out
    
    - *lppd*            bilevel lp
    - *GPRdict*         dictionary with (iso)enzymes mapped to reactions
    - *SubUdict*        dictionary with genes mapped to isoenzymes
    - *maxDelete*       maximum deletions allowed
    - *binGenes*        names of gene-binaries
    - *ReactionMap*     Enzyme-binaries mapped to Reaction-binaries
    - *EnzymeMap*       Gene-binaries mapped to Enzyme-binaries
    - *NoGene*          All binaries that correspond to multiple enzymes

    """
    binOn = []
    lplb = lppd.variables.get_lower_bounds()
    lpub = lppd.variables.get_upper_bounds()
    lpnam = lppd.variables.get_names()

    # Negative upperbounds and positive lower bounds cannot be in fused primal objective,
    # therefore, they are not included to be knocked out.
    for v_ in range(len(lpnam)):
        if lplb[v_] > 0:
            binOn.append(lpnam[v_])

        if lpub[v_] < 0:
            binOn.append(lpnam[v_])

    # Create binaries
    binGenes = []
    EnzymeMap = {}
    ReactionMap = {}
    for reaction in GPRdict:
        lppd.variables.add([0], [0], [1], types='B', names=['bin_'+reaction])
        iso = 0
        ReactionMap['bin_'+ reaction] = []
        for enzyme in GPRdict.get(reaction):
            iso += 1
            lppd.variables.add([0], [0], [1], types='B', names=['bin_E'+ str(iso) +'_'+reaction])
            ReactionMap.get('bin_' + reaction).append('bin_E'+str(iso) + '_' +reaction)

            EnzymeMap['bin_E'+str(iso) +'_'+reaction] = []
            for gene in enzyme:
                lppd.variables.add([0], [0], [1], types='B', names=['bin_'+gene])
                EnzymeMap.get('bin_E'+str(iso) +'_'+reaction).append('bin_'+gene)
                binGenes.append('bin_'+gene)

    NoGene = []
    for subunit in SubUdict:
        NoGene.append('bin_'+subunit)
        for gene in SubUdict.get(subunit):
            lppd.variables.add([0], [0], [1], types='B', names=['bin_'+gene])
            binGenes.append('bin_'+gene)

    binGenes = list(set(binGenes) - set(NoGene))


    # add number of minimal active fluxes or genes
    lppd.linear_constraints.add([cplex.SparsePair(binGenes, numpy.ones(len(binGenes)))],
                                 "G", rhs=[len(binGenes)-maxDelete], names=['MinActG'])

    return lppd, binGenes, ReactionMap, EnzymeMap, NoGene, binOn


def binFluxInd(lppd, fix_lower_bounds, fix_upper_bounds, GPRdict, primalVariables, USE_GENE, binPrimal, binOn):
    """
    Function that couples a flux value to a reaction binary
    Can set extra constraints that have not been included in the primal/dual (not recommended)
    
    
    - *lppd*                bilevel lp
    - *fix_lower_bounds*    lower bounds that are added to the problem 
                            after the bilevel lp is formed
    - *fix_upper_bounds*    lower bounds that are added to the problem 
                            after the bilevel lp is formed
    - *GPRdict*             dictionary with (iso)enzymes mapped to reactions
    - *primalVariables*     reactions
    - *USE_GENE*            Boolean indicator of use of gene-knockouts
    - *binPrimal*           Primal variables with binaries


    """

    # Reactions can only be switched off if they have genes
    if USE_GENE:
        for r_ in GPRdict:
            lppd.indicator_constraints.add([[r_], [1]], sense='E', rhs=0.0,
                                               indvar='bin_'+r_, complemented=1, name='ind_'+r_)
            binPrimal.append('bin_'+r_)

        # Reactions with a minimum flux have to stay on.
        for r_ in binOn:
            try:
                lppd.linear_constraints.add([cplex.SparsePair(['bin_'+r_], [1])], 'E', rhs=[1], names=['ON_bin_'+r_])
            except:
                print r_ +' does not have a binary'


    else:
        for r_ in primalVariables:
            if 'bin_'+ r_ in binPrimal:
                lppd.indicator_constraints.add([[r_], [1]], sense='E', rhs=0.0,
                                               indvar='bin_'+r_, complemented=1, name='ind_'+r_)

    return lppd, binPrimal

def dualConsInd(lppd, dualConstraints, binPrimal, dualVariables):
    """
        Function that constraints a dual constraint, if its associated primal variable has a binary
        
        - *lppd*            bilevel lp
        - *dualConstraints* names dual Constraints
        - *binPrimal*       primal variablies with a binary
        
        
        """
    ## adding dual binary constraints
    for c_ in dualConstraints:
        if 'bin_'+ c_ in binPrimal:
            c = 'd'+c_
            bVar = 'bin_'+c_
            spc = lppd.linear_constraints.get_rows(c)
            rhs = lppd.linear_constraints.get_rhs(c)
            spcX = [[lppd.variables.get_names(a) for a in spc.ind], spc.val]
            spcS = lppd.linear_constraints.get_senses(c)
            lppd.linear_constraints.delete(c)
            lppd.indicator_constraints.add(cplex.SparsePair(spcX[0], spcX[1]), sense=spcS, rhs=rhs, indvar=bVar, complemented=0, name='ind_'+c)
            names = []
            values = []
            for n_ in range(len(spc.ind)):
                n = lppd.variables.get_names(spc.ind[n_])
                if n not in dualVariables:
                    names.append(n)
                    values.append(spc.val[n_])
    return lppd

def andOrGPR(lppd, ReactionMap, EnzymeMap, SubUdict):
    """
    Function couples binaries of genes, (iso)enzymes and reactions, by defining constraints
    that define the AND and OR relationships.
    As a result, a reaction can be on, if one enzyme is present (OR). An enzyme can only be present
    if all the genes are on (AND).
    
    - *lppd*            bilevel lp
    - *ReactionMap*     Enzyme-binaries mapped to Reaction-binaries
    - *EnzymeMap*       Gene-binaries mapped to Enzyme-binaries
    - *SubUdict*        dictionary with genes mapped to isoenzymes


    """
    # add constrains to define the or/and relationships
    for reaction in ReactionMap:
        allVariables = []
        allVariables.append(reaction)
        for enzyme in ReactionMap.get(reaction):
            lppd.linear_constraints.add([cplex.SparsePair([reaction, enzyme], [1, -1])],
                                                 "G", rhs=[0], names=[enzyme])
            allVariables.append(enzyme)

        coefficients = [1] + [-1] *(len(allVariables)-1)
        lppd.linear_constraints.add([cplex.SparsePair(allVariables, coefficients)],
                                                 "L", rhs=[0], names=[reaction + '_E'])

    for enzyme in EnzymeMap:
        geneCount = 0
        allVariables = []
        allVariables.append(enzyme)
        for gene in EnzymeMap.get(enzyme):
            geneCount += 1
            lppd.linear_constraints.add([cplex.SparsePair([enzyme, gene], [1, -1])],
                                                 "L", rhs=[0], names=[gene + '_' + enzyme])
            allVariables.append(gene)

        coefficients = [1] + [-1]*(len(allVariables)-1)
        rest = geneCount*-1 + 1
        try:
            lppd.linear_constraints.add([cplex.SparsePair(allVariables, coefficients)],
                                                 "G", rhs=[rest], names=[gene + '_E'])
        except:
            print allVariables

    for subunit in SubUdict:
        allVariables = []
        allVariables.append('bin_'+subunit)
        for gene in SubUdict.get(subunit):
            lppd.linear_constraints.add([cplex.SparsePair(['bin_'+subunit, 'bin_'+ gene], [1, -1])],
                                                 "G", rhs=[0], names=[subunit+'_'+gene])
            allVariables.append('bin_'+gene)
        coefficients = [1] + [-1] *(len(allVariables)-1)
        lppd.linear_constraints.add([cplex.SparsePair(allVariables, coefficients)],
                                                 "L", rhs=[0], names=[subunit])

    return lppd


def useKnockWeight(lppd, binPrimal, bilevelObjective, binGenes, USE_GENE, KNOCKOUT_WEIGHTING_ALPHA):
    """
    Funtion adjusts the outer linear objective, by adding weight factors for all binaries
    E.g.: If the linear objective was R16, the new objective becomes R16 + 0.002 bin_G1 + 0.002 bin_G2 + 0.002 bin_G3 + ....
    In this way, the amount of binaries that are 1 are maximised, and thus the number of knockouts reduced. The weight factor 
    has a small value in order to keep R16 the ruling objective.
    
    - *lppd*            bilevel lp
    - *primalVariables  reactions
    - *bilevelObjective flux(es) of interest
    - *binGenes*        all gene binaries
    - *USE_GENE*        Boolean if genes are knocked out
    - *KNOOCKOUT_W_A*   Small value that is the weighting factor

    """

    weightObj = []
    KObinaries = []
    if USE_GENE:
        KObinaries = binGenes
    else:
        KObinaries = binPrimal

    weightObj.append(bilevelObjective)
    for b_ in KObinaries:
        # binaries are defined so 0 = off we want as many as possible on!
        # good because it makes the weighting less sensitive, would still keep alpha as small as possible
        weightObj.append((b_, KNOCKOUT_WEIGHTING_ALPHA))
    lppd.objective.set_linear(weightObj)

    del weightObj

    return lppd

def printLPsol(lp):
    """
    Prints the solution of an lp by variable in the terminal
    
    - *lp*              a cplex lp object
    """
    lp.solve()
    VariableSol = cbm.CBSolver.cplx_getOptimalSolution(lp)
    print('\nOptimum {}: {}'.format(VariableSol[1], VariableSol[2]))
    K = VariableSol[0].keys()
    K.sort()
    #for k_ in K:
    #   print(' {} : {}'.format(k_, VariableSol[0][k_]))

    return VariableSol[0]


def runOptKnock(modelFile, bilevelObjective, bio_reaction, objMinFactor, maxDelete, infinityValue, USE_GENE, USE_KNOCKOUT_WEIGHTING, KNOCKOUT_WEIGHTING_ALPHA, SOLUTION_FROM_OPTIMUM):


    time0 = time.time()

    # Routine for CPLEX method
  
    lp, lpName, GPRdict, SubUdict, Scale, maxDelete, infinityValue = readModel(modelFile, maxDelete, bio_reaction, infinityValue, USE_GENE)
    objVal, linCoeffVars = solveInner(lp, objMinFactor)
    lp, fix_lower_bounds, fix_upper_bounds = losenInnerObj(lp, infinityValue, linCoeffVars, objVal)
    #lp.write('02.lp')

    primalConstraints, primalVariables = characterizePrimal(lp)
    lpd = cplx_buildDual(lp)
    dualConstraints, dualVariables = characterizeDual(lpd,primalConstraints, lpName)

    lp.solve()
    lpd.solve()

    lppd = cplx_mergePrimalDual(lp, lpd, newobj=[bilevelObjective], osense='maximize')

    #lppd.write('04.lp')
    lppd.solve()
        
    if USE_GENE:
        binPrimal = []
        lppd, binGenes, ReactionMap, EnzymeMap, NoGene, binOn = createBinariesGenes(lppd, GPRdict, SubUdict, maxDelete)
    else:
        binOn = []
        binPrimal, lppd = createBinariesReac(lppd, primalVariables, maxDelete)
    #lppd.write('05.lp')

    lppd, binPrimal = binFluxInd(lppd, fix_lower_bounds, fix_upper_bounds, GPRdict, primalVariables, USE_GENE,  binPrimal, binOn)
    #lppd.write('06.lp')

    lppd = dualConsInd(lppd, dualConstraints, binPrimal, dualVariables)
    #lppd.write('07.lp')

    if USE_GENE:
        lppd = andOrGPR(lppd, ReactionMap, EnzymeMap, SubUdict)
    #lppd.write('08.lp')

    if USE_KNOCKOUT_WEIGHTING:
        if not USE_GENE:
            binGenes = []
        lppd = useKnockWeight(lppd, binPrimal, bilevelObjective, binGenes, USE_GENE, KNOCKOUT_WEIGHTING_ALPHA)


    # Write final MILP
    lppd.write(os.path.join('LPs', lpName+'_bilevel_milp.lp'))
    MILP = cplex.Cplex(os.path.join('LPs', lpName+'_bilevel_milp.lp'))

    # Get solution(s)
    MILP.parameters.mip.tolerances.mipgap.set(SOLUTION_FROM_OPTIMUM)

    try:
        MILP.solve()
        val1 = MILP.solution.get_objective_value()/Scale
        val2 = MILP.solution.get_values(bio_reaction)/Scale
        sol = printLPsol(MILP)
    except cplex.exceptions.CplexSolverError:
        print('\nFusedsolveFailure\n')

    time1 = time.time()

    delGen = []
    delReact = []

    if USE_GENE:
        if val1 != 0.0:
            delG = [x for x in sol if x.startswith('bin_GENE_') and sol[x] == 0.0]
            delG = set(delG) - set(NoGene)
            for g_ in delG:
                delGen.append(g_.replace('bin_GENE_',''))
    else:
        if val1 != 0.0:
            delReac = [x for x in sol if x.startswith('bin_R') and sol[x] == 0.0]
            for r_ in delReac:
                delReact.append(r_.replace('bin_',''))



    # Print results
    print('\npyOpt Solution\n=================')
    print('\nProduct reaction ({}): {}'.format(bilevelObjective[0], round(val1,5)))
    print('\nBiomass reaction ({}): {}'.format(bio_reaction, round(val2,5)))
    if USE_GENE:
        print('\nDeleted genes ({}): {}'.format(len(delGen), delGen))
    else:
        print('\nDeleted reactions ({}): {}'.format(len(delReact), delReact))
    print('\nMaximum deletions: {}'.format(maxDelete))
    print('\nUse knockout weighting (minimization): {}'.format(USE_KNOCKOUT_WEIGHTING))
    if USE_KNOCKOUT_WEIGHTING:
        print('WARNING! The outer optimum value is slightly lower, due to use of knockout weighting')
    print('\n=================\n')

    
    return [bilevelObjective, val1, bio_reaction, val2, delGen,delReact]
