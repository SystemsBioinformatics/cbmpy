import cbmpy

ucmod = cbmpy.readSBML3FBC('test_user_def_rc.xml')
ucmod.importUserConstraints('test_user_def_rc.xml.user_constraints.json')
print(ucmod.user_constraints)

UC = {'cons1': {'fluxes': [[1, 'R16'], [1, 'R15']], 'operator': 'LE', 'rhs': 1.0},
 'cons2': {'fluxes': [[1, 'R15']], 'operator': 'GE', 'rhs': -1},
 'cons3': {'fluxes': [[1, 'R14']], 'operator': 'E', 'rhs': 0.0}}

for u in UC:
    if UC[u]['operator'] == 'E':
        lb = UC[u]['rhs']
        ub = UC[u]['rhs']
    elif UC[u]['operator'] in ['G', 'GE']:
        lb = UC[u]['rhs']
        ub = cbmpy.INF
    elif UC[u]['operator'] in ['L', 'LE']:
        lb = cbmpy.NINF
        ub = UC[u]['rhs']

    components = []
    for c in UC[u]['fluxes']:
        components.append((c[0], c[1], 'linear'))

    print(lb, ub)

    ucnew = ucmod.createUserDefinedConstraint(u, lb, ub, components)
    ucmod.addUserDefinedConstraint(ucnew)
