import cplex
import numpy
import cbmpy
import time
import cbmpy.CBCPLEX as wcplex
from packaging import version as pkgver

from DCFBA.ToyModels.model_c import build_model_C

_Stime = time.time()

HAVE_SYMPY = False
try:
    import sympy

    if pkgver.parse(sympy.__version__) >= pkgver.Version("0.7.4"):
        HAVE_SYMPY = True
    else:
        del sympy
        print(
            "\nWARNING: SymPy version 0.7.5 or newer is required for symbolic matrix support."
        )
except ImportError:
    HAVE_SYMPY = False

HAVE_SCIPY = False
try:
    from scipy.sparse import csr

    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

CPLX_LP_PARAMETERS = {
    "simplex.tolerances.optimality": 1e-6,
    "simplex.tolerances.feasibility": 1e-6,
}


def cplex_constructProbfromFBA(fba, fname=None):
    # define model and add variables
    prob = cplex.Cplex()
    # define simplex tolerances for the model
    prob.parameters.simplex.tolerances.optimality.set(
        CPLX_LP_PARAMETERS["simplex.tolerances.optimality"]
    )
    prob.parameters.simplex.tolerances.feasibility.set(
        CPLX_LP_PARAMETERS["simplex.tolerances.feasibility"]
    )
    # print(lp.parameters.simplex.get_changed())
    prob.set_problem_name("%s" % (fba.getId()))
    prob.variables.add(
        names=fba.N.col,
        types=[prob.variables.type.continuous] * len(fba.N.col),
    )

    try:
        # define objective
        osense = fba.getActiveObjective().operation.lower()
        if osense in ["minimize", "minimise", "min"]:
            prob.objective.set_sense(prob.objective.sense.minimize)
        elif osense in ["maximize", "maximise", "max"]:
            prob.objective.set_sense(prob.objective.sense.maximize)
        else:
            raise RuntimeError(
                "\n%s - is not a valid objective operation" % osense
            )
        prob.objective.set_name(fba.getActiveObjective().getId())
    except AttributeError:
        print("\nWARNING(CPLEX create LP): no objective function defined")
    # TODO in the model active objective set a variable that says if
    # objective is quadratic !implementation needed!

    cplex_buildLinearConstraints(prob, fba, fname)

    cplex_constructLPfromFBA(prob, fba)
    cplex_constructQPfromFBA(prob, fba)

    return prob


def cplex_constructQPfromFBA(prob, fba):
    # target = prob.parameters.optimalitytarget.values
    # prob.parameters.optimalitytarget.set(target.optimal_global)

    cmat = []
    for fo in fba.getActiveObjective().QPObjective:
        v = fo[1]
        if fo[0][0] == fo[0][1]:
            v *= 2

        cmat.append((fo[0][0], fo[0][1], v))

    prob.objective.set_quadratic_coefficients(cmat)

    return prob


def cplex_constructLPfromFBA(prob, fba):
    prob.objective.set_linear(
        [
            (fo.reaction, fo.coefficient)
            for fo in fba.getActiveObjective().fluxObjectives
        ]
    )


def cplex_buildLinearConstraints(prob, fba, fname):
    lin_expr = []
    rhs = []
    names = []
    senses = []
    lp = prob
    if HAVE_SYMPY and fba.N.__array_type__ == sympy.MutableDenseMatrix:
        print("INFO: CPLEX requires floating point, converting N")
        Nmat = numpy.array(fba.N.array).astype("float")
        RHSmat = numpy.array(fba.N.RHS).astype("float")
        if fba.CM != None:
            CMmat = numpy.array(fba.CM.array).astype("float")
            CMrhs = numpy.array(fba.CM.RHS).astype("float")
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
        lin_expr.append(
            cplex.SparsePair(
                [fba.N.col[c] for c in nz], [Nmat[n, c] for c in nz]
            )
        )
        rhs.append(RHSmat[n])
        senses.append(wcplex.cplx_fixConSense(fba.N.operators[n]))
        names.append(fba.N.row[n])
    # print senses
    lp.linear_constraints.add(
        lin_expr=lin_expr, senses=senses, rhs=rhs, names=names
    )
    # print 'New style lc:', time.time() - tnew

    # add user defined constraints
    if fba.CM != None:
        # the numpy way
        # t2new = time.time()
        lin_expr = []
        rhs = []
        names = []
        senses = []
        for n in range(CMmat.shape[0]):
            if not GOSCI:
                nz = CMmat[n].nonzero()[0]
            else:
                nz = CMmat[n].nonzero()[1]
            lin_expr.append(
                cplex.SparsePair(
                    [fba.CM.col[c] for c in nz], [CMmat[n, c] for c in nz]
                )
            )
            rhs.append(CMrhs[n])
            senses.append(wcplex.cplx_fixConSense(fba.CM.operators[n]))
            names.append(fba.CM.row[n])
        lp.linear_constraints.add(
            lin_expr=lin_expr, senses=senses, rhs=rhs, names=names
        )
        # print 'New style lc:', time.time() - t2new

    # add bounds
    lb = []
    ub = []
    for b_ in fba.flux_bounds:
        btype = b_.getType()
        bvalue = b_.getValue()
        if bvalue in ["Infinity", "inf", "Inf", "infinity"]:
            bvalue = cplex.infinity
        elif bvalue in ["-Infinity", "-inf", "-Inf", "-infinity"]:
            bvalue = -cplex.infinity
        elif numpy.isinf(bvalue):
            if bvalue > 0.0:
                bvalue = cplex.infinity
            else:
                bvalue = -cplex.infinity
        if btype == "lower":
            lb.append((b_.reaction, bvalue))
        elif btype == "upper":
            ub.append((b_.reaction, bvalue))
        elif btype == "equality":
            ub.append((b_.reaction, bvalue))
            lb.append((b_.reaction, bvalue))
    if len(lb) > 0:
        lp.variables.set_lower_bounds(lb)
    if len(ub) > 0:
        lp.variables.set_upper_bounds(ub)

    print("\ncplx_constructLPfromFBA time: {}\n".format(time.time() - _Stime))
    if fname != None:
        lp.write(fname + ".lp", filetype="lp")
    return lp


# m_c = build_model_C()
# m_c.getReaction("S_exchange").setLowerBound(-100)

# obj = m_c.getActiveObjective()
# obj.QPObjective = [(["BM_e_C_exchange", "BM_e_C_exchange"], 1.0)]

# m_c.buildStoichMatrix()
# cbmpy.doFBA(m_c)
# cplex_constructPorbfromFBA(m_c)
