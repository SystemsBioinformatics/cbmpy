import cbmpy

# read my first SBML3 FBCv3 file

print('\nREADING ...\n')
cmod3 = cbmpy.CBRead.readSBML3FBC('the_third.fbc3.xml')


print('FOdata:', cmod3.getActiveObjective().getFluxObjectiveData())
print('Linear objective variables:', cmod3.getActiveObjective().getLinearVariables())
print('Quadratic objective variables:', cmod3.getActiveObjective().getQuadraticVariables()[0].getReactionId())

print('\nWRITING ...\n')
cbmpy.CBWrite.writeSBML3FBCV3(cmod3, 'the_virtually_finished.fbc3.xml', compress_bounds=True)

print('done ....')
