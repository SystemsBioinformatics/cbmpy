import cbmpy

# read my first SBML3 FBCv3 file

print('\nREADING ...\n')
cmod3 = cbmpy.CBRead.readSBML3FBC('the_third.fbc3.xml')


print('\nWRITING ...\n')
cbmpy.CBWrite.writeSBML3FBCV3(cmod3, 'the_fourth-firstroundtrip.fbc3.xml', compress_bounds=True)

print('done.')