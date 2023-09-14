import cbmpy

# read my first SBML3 FBCv3 file

print('\nOpen the_xenomorph.fbc3.xml to see how objectives and constraints are encoded.\n')
# cmod3 = cbmpy.CBRead.readSBML3FBC('test_user_def_rc.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('e_coli_core.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('cbmpy_test_ecoli')

cmod3 = cbmpy.CBRead.readSBML3FBC('the_xenomorph.fbc3.xml')

# Some methods of interest, there are more than these so explore them a bit

# constraints
cmod3.createUserDefinedConstraint
cmod3.user_defined_constraints
cbmpy.CBModel.UserDefinedConstraint.createConstraintComponent
cbmpy.CBModel.UserDefinedConstraint.getConstraintComponent
cbmpy.CBModel.UserDefinedConstraint.getConstraintComponentForVariable
# convert old style constraints to new
cmod3.convertUserConstraintsToUserDefinedConstraints

# objectives
cmod3.getActiveObjective().getLinearFluxObjectives
cmod3.getActiveObjective().getQuadraticFluxObjectives
cmod3.getActiveObjective().getQuadraticBivariateFluxObjectives


cbmpy.doFBA(cmod3)

print('Questions\nCan you encode a QP model that actually works, do I need to change any of the CPLEX options for a QP?')


# # this is all to test bounds, was used to debug stuff

# cmod3.createReaction("R_VN", "Import S external to the cytosol", reversible=True)
# cmod3.createReactionReagent("R_VN", cmod3.species[3].getId(), -1)
# cmod3.createReactionReagent("R_VN", cmod3.species[4].getId(), 1)
# cmod3.setReactionBounds("R_VN", -1000.0, 1000.0)
# print(cmod3.getReactionBounds("R_VN"))