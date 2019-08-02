import os, time, numpy
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import cbmpy as cbm

print(os.sys.argv)

mBase = os.sys.argv[2]
testDir = os.path.join(os.sys.argv[1], mBase)
outDir = os.sys.argv[3]

if not os.path.exists(outDir):
    os.makedirs(outDir)

cmod = cbm.CBRead.readSBML3FBC('%s-sbml-l3v1.xml' % mBase, testDir)
## cbm.CBSolver.cplx_analyzeModel(cmod, lpFname=os.path.join(outDir, mBase))
## cbm.CBSolver.glpk_analyzeModel(cmod)
cbm.CBSolver.analyzeModel(cmod)

F = open(os.path.join(testDir, '%s-settings.txt' % mBase), 'r')
variables = None
for l_ in F:
    if l_.startswith('variables:'):
        variables = l_.strip().split(':')[1].strip().split(',')
        break
#print variables
data = []
for v_ in variables:
    vtest = cmod.getReaction(v_)
    if vtest == None:
        data.append(cmod.getActiveObjective().getValue())
    else:
        data.append(cmod.getReaction(v_).getValue())
#print data

F = open(os.path.join(outDir, mBase)+'.csv', 'w')
out = ''
for v_ in variables:
    out += '%s,' % v_
out = '%s\n' % out[:-1]
for d_ in data:
    out += '%s,' % d_
out = '%s\n' % out[:-1]
F.write(out)
F.close()

