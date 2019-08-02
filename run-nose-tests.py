import os, zipfile, cbmpy, libsbml
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import nose

__test_dir__ = os.path.join(cDir, 'nosetests')
test_set = 1
test_bundle = os.path.join(__test_dir__, 'cbmpy_test_bundle_v{}.zip.py'.format(test_set))
test_file = os.path.join(__test_dir__, 'installed.v{}'.format(test_set))
if not os.path.exists(test_file):
    print('Installing test files (v{})...'.format(test_set))
    zfile = zipfile.ZipFile(test_bundle, allowZip64=True)
    zfile.extractall(path=__test_dir__)
    zfile.close()
    #os.remove(test_bundle)
    del zipfile, zfile

os.chdir(__test_dir__)

nose.run()
