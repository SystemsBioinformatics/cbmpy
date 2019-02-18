try:
    __test_dir__ = os.path.join(__CBCONFIG__['CBMPY_DIR'], 'nosetests')
    test_set = 1
    test_bundle = os.path.join(__test_dir__, 'cbmpy_test_bundle_v{}.zip.py'.format(test_set))
    test_file = os.path.join(__test_dir__, 'installed.v{}'.format(test_set))
    if not os.path.exists(test_file):
        import zipfile
        print('Installing test files (v{})...'.format(test_set))
        zfile = zipfile.ZipFile(test_bundle, allowZip64=True)
        zfile.extractall(path=__test_dir__)
        zfile.close()
        #os.remove(test_bundle)
        del zipfile, zfile
    from . import nosetests
    test = nosetests.run
    del nosetests, test_bundle, test_file, test_set
except ImportError:
    def test():
        print("ERROR: Tests not installed")