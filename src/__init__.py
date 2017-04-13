"""
CBMPy: __init__.py
==================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2017 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: __init__.py 575 2017-04-13 12:18:44Z bgoli $)

"""
#
# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os

# Setup some environment
__SILENT_START__ = False
from . import CBConfig
__CBCONFIG__ = CBConfig.__CBCONFIG__
__version__ = CBConfig.__CBCONFIG__['VERSION']
__CBCONFIG__['CBMPY_DIR'] = os.path.split(CBConfig.__file__)[0]

# This is just a hack for backwards compatibility with existing scripts
# that enables "from cbmpy import *" to work - bgoli

# lets add some shortcuts that will be core commands
__HAVE_SOLVER__ = False
try:
    from . import CBSolver
    analyzeModel = CBSolver.analyzeModel
    FluxVariabilityAnalysis = CBSolver.FluxVariabilityAnalysis
    MinimizeSumOfAbsFluxes = CBSolver.MinimizeSumOfAbsFluxes
    doFBA = CBSolver.analyzeModel
    doFVA = CBSolver.FluxVariabilityAnalysis
    doFBAMinSum = CBSolver.MinimizeSumOfAbsFluxes
    __HAVE_SOLVER__ = True
except (ImportError, AttributeError):
    print('No solver present, unable to create shortcuts')

from .CBWrite import writeFVAtoCSV, writeModelToExcel97
from .CBXML import _HAVE_SBML_
if _HAVE_SBML_:
    from .CBXML import sbml_fileFindVersion, sbml_fileValidate
    from .CBRead import readSBML3FBC, readSBML2FBA, readCOBRASBML
    from .CBWrite import writeSBML3FBC, writeModelToExcel97, writeModelToCOMBINEarchive, writeCOBRASBML, writeSBML3FBCV2
else:
    print('\nWARNING: No SBML support, top-level SBML read/write functions disabled.\n')

# CBMPy GUI development has ended with WxGUI, download Metatoolkit instead
try:
    from .CBWx import runModelEditor as loadCBGUI
except ImportError as ex:
    print('GUI module not available')

##import CBMultiEnv to be deprecated
#try:
    #from . import CBMultiEnv
#except ImportError as ex:
    #print(ex)
    #print('MultiEnvironment module not available')
    #CBMultiEnv = None
#import CBMultiCore
try:
    from . import CBMultiCore
except ImportError as ex:
    print(ex)
    print('MultiCore module not available')
    CBMultiCore = None

INF = float('inf')
NINF = -float('inf')
NAN = float('nan')

from . import CBModel, CBDataStruct, CBModelTools, CBRead, CBReadtxt, CBTools, CBWrite, CBXML, CBNetDB, CBPlot, CBMultiModel
from .CBTools import deSerialize, deSerializeFromDisk
from . import PyscesSED
SED = PyscesSED.SEDCBMPY

try:
    from . import nosetests
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
    test = nosetests.run
    del nosetests, test_bundle, test_file, test_set
except ImportError:
    def test():
        print("ERROR: Tests not installed")

del os

if not __SILENT_START__:
    print('\nCBMPy environment\n******************')
    rev = CBConfig.__CBCONFIG__['VERSION_STATUS']
    #print('Release: {}'.format(__version__.replace('.{}'.format(rev),'')))
    print('Revision: r{}'.format(CBConfig.__CBCONFIG__['VERSION_STATUS']))
    #print 'pysces.model_dir = ' + model_dir
    #print 'pysces.output_dir = ' + output_dir

    print('\n\n***********************************************************************')
    print('* Welcome to CBMPy ({}) - PySCeS Constraint Based Modelling        *'.format(__version__.replace('.{}'.format(rev),'')))
    print('*                http://cbmpy.sourceforge.net                         *')
    print('* Copyright(C) Brett G. Olivier 2014 - 2017                           *')
    print('* Dept. of Systems Bioinformatics                                     *')
    print('* Vrije Universiteit Amsterdam, Amsterdam, The Netherlands            *')
    print('* CBMPy is developed as part of the BeBasic MetaToolKit Project       *')
    #print('* PySCeS is distributed under the PySCeS (BSD style) licence, see     *')
    print('* Distributed under the GNU GPL v 3.0 licence, see                    *')
    print('* LICENCE (supplied with this release) for details                    *')
    ##  print('*                 ** Read about PySCeS **                             *')
    ##  print('* http://bioinformatics.oupjournals.org/cgi/content/short/21/4/560    *')
    print('***********************************************************************\n')

