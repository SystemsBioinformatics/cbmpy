from cobra.io import load_model
import cbmpy
from cbmpy.CBModel import Model, Objective
from testingQP import cplex_constructProbfromFBA
import cbmpy.CBCPLEX as wcplex

model = load_model("textbook")
model.solver = "cplex"
sum_two = model.problem.Constraint(
    model.reactions.FBA.flux_expression + model.reactions.NH4t.flux_expression,
    lb=2,
    ub=2,
)
model.add_cons_vars(sum_two)

quadratic_objective = model.problem.Objective(
    0.5 * model.reactions.NH4t.flux_expression**2
    + 0.5 * model.reactions.FBA.flux_expression**2
    - model.reactions.FBA.flux_expression,
    direction="min",
)
model.objective = quadratic_objective
solution = model.optimize(objective_sense=None)

print(solution.fluxes["NH4t"], solution.fluxes["FBA"])


print(solution)

# cbmpy recreate:
model: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")


model.addUserConstraint("sum2", [[1, "R_FBA"], [1, "R_NH4t"]], "=", 2)
model.createObjectiveFunction("R_FBA", -1, "minimize")

obj: Objective = model.getActiveObjective()

obj.QPObjective = []
obj.QPObjective.append((("R_NH4t", "R_NH4t"), 0.5))
obj.QPObjective.append((("R_FBA", "R_FBA"), 0.5))

model.buildStoichMatrix()

prob = cplex_constructProbfromFBA(model)
print(" quad  = ", prob.objective.get_quadratic())
print(" lin   = ", prob.objective.get_linear())
print(" sense = ", prob.objective.get_sense())


prob.solve()
print(prob.solution.status[prob.solution.get_status()])
print("Solution value  = ", prob.solution.get_objective_value())
sol, objname, objval = wcplex.cplx_getOptimalSolution(prob)

(
    model.objectives[model.activeObjIdx].solution,
    model.objectives[model.activeObjIdx].value,
) = (sol, objval)
for r in model.reactions:
    rid = r.getId()
    if rid in sol:
        r.value = sol[rid]
    else:
        r.value = None

FBAsol = model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
print(FBAsol["R_FBA"])
print(FBAsol["R_NH4t"])
