import cbmpy

ucmod = cbmpy.readSBML3FBC('test_user_def_rc.xml')
ucmod.importUserConstraints('test_user_def_rc.xml.user_constraints.json')
print(ucmod.user_constraints)

ucmod.convertUserConstraintsToUserDefinedConstraints()
print(ucmod.user_constraints, ucmod.user_defined_constraints)
print(list(ucmod.__global_id__.keys()))