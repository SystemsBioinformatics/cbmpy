import os, time, numpy
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pyscescbm as cbm

model_name = '130517BP_NoTCA.xml'
cmod = cbm.CBRead.readSBML2FBA(os.path.join(cDir, 'models', model_name))
cmodS = cmod.clone()

cbm.analyzeModel(cmod)
fva, fvan = cbm.CBMultiCore.runMultiCoreFVA(cmod, procs=4)
cbm.CBWrite.writeFVAtoCSV(model_name.replace('.xml',''), fva, fvan, fbaObj=cmod)
cbm.CBWrite.writeModelToExcel97(cmod, model_name.replace('.xml',''))

exch_reactions = cmod.getReactionIds('R_EX_')
dump = []
# reactions are defined back to front so fwd is "o" and and backward is "i"
for xr_ in exch_reactions:
    dump.append(cbm.CBTools.splitSingleReversibleReaction(cmodS, xr_, xr_.replace('R_EX_', 'R_EXo_'),\
                                                          xr_.replace('R_EX_', 'R_EXi_')))
cbm.analyzeModel(cmodS)
fva, fvan = cbm.CBMultiCore.runMultiCoreFVA(cmodS, procs=4)
cbm.CBWrite.writeFVAtoCSV(model_name.replace('.xml','.split'), fva, fvan, fbaObj=cmodS)
cbm.CBWrite.writeSBML2FBA(cmodS, model_name.replace('.xml','.split.xml'))
cbm.CBWrite.writeModelToExcel97(cmodS, model_name.replace('.xml','.split'))

print('\nObjectives')
print('base model: {}'.format(cmod.getObjFuncValue()))
print('split model: {}'.format(cmodS.getObjFuncValue()))