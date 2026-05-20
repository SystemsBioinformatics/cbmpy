import cbmpy

# read my first SBML3 FBCv3 file

print('\nOpen the_xenomorph.fbc3.xml to see how objectives and constraints are encoded.')
print('Open combined_e_coli_core.fbc3.xml to see how quadratric objectives are encoded and solved.\n')
# cmod3 = cbmpy.CBRead.readSBML3FBC('test_user_def_rc.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('e_coli_core.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('cbmpy_test_ecoli')
# cmod3 = cbmpy.CBRead.readSBML3FBC('the_xenomorph.fbc3.xml')
# cmod3 = cbmpy.readSBML3FBC('test_community.fbc3.xml')
cmod3 = cbmpy.readSBML3FBC('combined_e_coli_core.fbc3.xml')

# # Some methods of interest, there are more than these so explore them a bit

# # constraints
# cmod3.createUserDefinedConstraint
# cmod3.user_defined_constraints
# cbmpy.CBModel.UserDefinedConstraint.createConstraintComponent
# cbmpy.CBModel.UserDefinedConstraint.getConstraintComponent
# cbmpy.CBModel.UserDefinedConstraint.getConstraintComponentForVariable
# # convert old style constraints to new
# cmod3.convertUserConstraintsToUserDefinedConstraints

# objective information
print('\nobjective information')
print('FOdata:', cmod3.getActiveObjective().getFluxObjectiveData())
print('Linear objective variables:', cmod3.getActiveObjective().getLinearFluxObjectives())
print('Quadratic objective variables:', cmod3.getActiveObjective().getQuadraticFluxObjectives())
print('Quadratic Bivariate objective variables:', cmod3.getActiveObjective().getQuadraticBivariateFluxObjectives())


cbmpy.doFBA(cmod3)
# cpx = cbmpy.CBCPLEX.cplx_analyzeModel(cmod3, return_lp_obj=True, lpFname='combined_e_coli_core.fbc3')

