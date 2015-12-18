import os, time, numpy, subprocess, zipfile
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pyscescbm as cbm


testDir = os.path.join(cDir, 'tests')
if not os.path.exists(testDir):
    os.makedirs(testDir)
    if os.path.exists(os.path.join(cDir, 'fbc_test_suite.zip')):
        print('Unzipping Test suite to: '.format(testDir))
        zip_ref = zipfile.ZipFile(os.path.join(cDir, 'fbc_test_suite.zip'), 'r')
        zip_ref.extractall(testDir)
        zip_ref.close()


outDir = os.path.join(cDir, 'fbc_test_suite_out')
if not os.path.exists(outDir):
    os.makedirs(outDir)

tests = os.listdir(testDir)
print tests

for T in tests:
    retcode = subprocess.call([os.sys.executable, 'run_fbc_test.py', testDir, T, outDir])

testResults = []
for T in tests:
    rF = file(os.path.join(testDir, T, '%s-results.csv' % T), 'r')
    r0 = rF.read().strip()
    r0 = [a.split(',') for a in r0.split('\n')]
    r00 = {}
    for j_ in range(len(r0[0])):
        r00[r0[0][j_]] = r0[1][j_]
    rF.close()
    rF2 = file(os.path.join(outDir, '%s.csv' % T), 'r')
    r1 = rF2.read().strip()
    r1 = [a.split(',') for a in r1.split('\n')]
    r11 = {}
    for j_ in range(len(r1[0])):
        r11[r1[0][j_]] = r1[1][j_]
    rF2.close()
    testRes = 'fail'
    for k_ in r00:
        if k_ in r11:
            if r00[k_].upper() == r11[k_].upper():
                testRes = 'pass'
    testResults.append((T, testRes))

print '\n\nFBC Test Suite Results:\n===================='
for t_ in testResults:
    print 'Test %s: %s' % (t_[0], t_[1])
