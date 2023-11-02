import os, pprint

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import cbmpy

fbase = cbmpy.readSBML3FBC('cbmpy_test_core')

# DEFINE TEST MODEL
test = fbase.clone()
test.addUserConstraint('cons1', fluxes=[(1, 'R16'), (1, 'R15')], operator='=', rhs=0.8)
test.addUserConstraint('cons2', fluxes=[(1, 'R15')], operator='>=', rhs=0.1)
test.addUserConstraint('cons3', fluxes=[(1, 'R14')], operator='>=', rhs=0.1)
test.addUserConstraint('cons4', fluxes=[(1, 'R15')], operator='<=', rhs=1.0)

lp = cbmpy.CBSolver.cplx_analyzeModel(
    test,
    lpFname='dbg',
    return_lp_obj=True,
    with_reduced_costs='scaled',
    with_sensitivity=False,
    del_intermediate=False,
    build_n=True,
    quiet=False,
    oldlpgen=False,
    method='o',
)

OBJSENSE, RHSSENSE, BNDSENSE = cbmpy.CBSolver.cplx_getSensitivities(lp)

print('\nObjective sensitivity (reduced cost)')
pprint.pprint(OBJSENSE)
print('\nRight hand side sensitivity (might be related to a shadow price)')
pprint.pprint(RHSSENSE)
print('\nBound sensitivity')
pprint.pprint(BNDSENSE)

cbmpy.writeSBML3FBCV2(test, 'test_user_def_rc.xml')