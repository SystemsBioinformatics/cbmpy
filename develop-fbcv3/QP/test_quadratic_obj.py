import cbmpy

# read my first SBML3 FBCv3 file

print('\nREADING ...\n')
cmod3 = cbmpy.CBRead.readSBML3FBC('e_coli_core.xml')
# cmod3 = cbmpy.CBRead.readSBML3FBC('test_user_def_rc.xml')

rid =  cmod3.reactions[20].getId()

print('getReactionBounds')
print('Known ({}):'.format(rid), cmod3.getReactionBounds(rid))
print('UnKnown:', cmod3.getReactionBounds('Fred'))

print('getUpperBound')
print('Known ({}):'.format(rid), cmod3.getReaction(rid).getLowerBound())
print('Known ({}):'.format(rid), cmod3.getReaction(rid).getUpperBound())
try:
    print('UnKnown (this should fail):', cmod3.getReaction('fred').getUpperBound())
except AttributeError:
    print('UnKnown (this should fail)')
