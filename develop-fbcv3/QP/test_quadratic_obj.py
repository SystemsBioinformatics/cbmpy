import cbmpy

# read my first SBML3 FBCv3 file

print('\nREADING ...\n')
# cmod3 = cbmpy.CBRead.readSBML3FBC('test_user_def_rc.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('e_coli_core.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('cbmpy_test_ecoli')

cmod3 = cbmpy.CBRead.readSBML3FBC('the_xenomorph.fbc3.xml')

# cmod3.getActiveObjective().setOperation('minimize')

cbmpy.doFBA(cmod3)


# # this is all to test bounds, was used to fix Stevens bugs
# rid =  cmod3.reactions[20].getId()
# print('getReactionBounds')
# print('Known ({}):'.format(rid), cmod3.getReactionBounds(rid))
# print('UnKnown:', cmod3.getReactionBounds('Fred'))

# print('getUpperBound')
# print('Known ({}):'.format(rid), cmod3.getReaction(rid).getLowerBound())
# print('Known ({}):'.format(rid), cmod3.getReaction(rid).getUpperBound())
# try:
    # print('UnKnown (this should fail):', cmod3.getReaction('fred').getUpperBound())
# except AttributeError:
    # print('UnKnown (this should fail)')


# R = cbmpy.CBModel.Reaction('newR')
# cmod3.addReaction(R, create_default_bounds=True, silent=True)

# for r in cmod3.reactions:
    # lb1 = r.getLowerBound()
    # ub1 = r.getUpperBound()
    # b = cmod3.getReactionBounds(r.getId())
    # lb2 = cmod3.getReactionLowerBound(r.getId())
    # ub2 = cmod3.getReactionUpperBound(r.getId())
    # print(r.getId(), b)
    # print('lower pass:', lb1 == lb2)
    # print('upper pass:', ub1 == ub2)

# cmod3.createReaction("R_VN", "Import S external to the cytosol", reversible=True)
# cmod3.createReactionReagent("R_VN", cmod3.species[3].getId(), -1)
# cmod3.createReactionReagent("R_VN", cmod3.species[4].getId(), 1)
# cmod3.setReactionBounds("R_VN", -1000.0, 1000.0)
# print(cmod3.getReactionBounds("R_VN"))