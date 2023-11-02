import cbmpy

# read my first SBML3 FBCv3 file

print('\nREADING ...\n')
cmod3 = cbmpy.readSBML3FBC('the_third.fbc3.xml')


print('FOdata:', cmod3.getActiveObjective().getFluxObjectiveData())
print('Linear objective variables:', cmod3.getActiveObjective().getLinearFluxObjectives())
print('Quadratic objective variables:', cmod3.getActiveObjective().getQuadraticFluxObjectives())
print('Quadratic Bivariate objective variables:', cmod3.getActiveObjective().getQuadraticBivariateFluxObjectives())

print('\nWRITING ...\n')
cbmpy.writeSBML3FBCV3(cmod3, 'the_virtually_finished.fbc3.xml', compress_bounds=True)

print('done ....')
