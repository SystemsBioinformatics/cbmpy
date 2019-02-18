# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, nose
import cbmpy

fDir = os.path.split(__file__)[0]
#fDir = os.path.split(fDir)[0]
print('TestBase: ', fDir)

def run():
    cDir = os.getcwd()
    os.chdir(fDir)
    res = nose.run()
    os.chdir(cDir)
    if res:
        print('\nAll tests passed.\n')
    else:
        print('\nPlease ignore the second test failure:\n\n\"ERROR: Failure: ImportError (No module named cbmpy)\"\n\nthis is a known glitch in the testing framework that is currently being worked on.\n')

class TestBase:
    def test_import(self):
        import cbmpy
