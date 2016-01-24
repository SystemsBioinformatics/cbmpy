import os, nose
fDir = os.path.split(__file__)[0]

def run():
    cDir = os.getcwd()
    os.chdir(fDir)
    res = nose.run()
    os.chdir(cDir)
    if res:
        print('\nAll tests passed.\n')
