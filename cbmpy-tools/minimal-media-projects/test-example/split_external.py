import os, time, numpy
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import cbmpy as cbm

model_name = 'iAF692.xml'
cmod = cbm.readSBML3FBC(os.path.join(cDir, 'models', model_name))
cmodS = cmod.clone()

cbm.doFBA(cmod)
fva, fvan = cbm.CBMultiCore.runMultiCoreFVA(cmod, procs=4)
cbm.CBWrite.writeFVAtoCSV(fva, fvan, model_name.replace('.xml',''), fbaObj=cmod)
cbm.writeModelToExcel97(cmod, model_name.replace('.xml',''))

exch_reactions = cmod.getReactionIds('R_EX_')
dump = []
# reactions are defined back to front so fwd is "o" and and backward is "i"
for xr_ in exch_reactions:
    dump.append(cbm.CBTools.splitSingleReversibleReaction(cmodS, xr_, xr_.replace('R_EX_', 'R_EXo_'),\
                                                          xr_.replace('R_EX_', 'R_EXi_')))
cbm.doFBA(cmodS)
fva, fvan = cbm.CBMultiCore.runMultiCoreFVA(cmodS, procs=4)
cbm.CBWrite.writeFVAtoCSV(fva, fvan, model_name.replace('.xml','.split'), fbaObj=cmodS)
cbm.writeSBML3FBC(cmodS, model_name.replace('.xml','.split.xml'))
cbm.CBWrite.writeModelToExcel97(cmodS, model_name.replace('.xml','.split'))

print('\nObjectives')
print('base model: {}'.format(cmod.getObjFuncValue()))
print('split model: {}'.format(cmodS.getObjFuncValue()))