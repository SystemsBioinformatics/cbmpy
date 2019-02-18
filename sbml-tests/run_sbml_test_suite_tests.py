import os, subprocess
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

os.chdir(os.path.join(cDir, 'sbml_test_suite'))
retcode = subprocess.call([os.sys.executable, 'run_all_tests.py'])

os.chdir(cDir)
os.sys.exit(retcode)
