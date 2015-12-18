"""
CBMPy: __init__.py
==================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2015 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: __init__.py 358 2015-08-11 10:41:01Z bgoli $)

"""
##
# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

# Setup some environment
__SILENT_START__ = False
from . import CBConfig
__CBCONFIG__ = CBConfig.__CBCONFIG__
__version__ = CBConfig.__CBCONFIG__['VERSION']


# This is just a hack for backwards compatibility with existing scripts
# that enables "from pyscescbm import *" to work - bgoli

# lets add some shortcuts that will be core commands
try:
    from . import CBSolver
    analyzeModel = CBSolver.analyzeModel
    FluxVariabilityAnalysis = CBSolver.FluxVariabilityAnalysis
    MinimizeSumOfAbsFluxes = CBSolver.MinimizeSumOfAbsFluxes
    doFBA = CBSolver.analyzeModel
    doFVA = CBSolver.FluxVariabilityAnalysis
    doFBAMinSum = CBSolver.MinimizeSumOfAbsFluxes
except ImportError:
    print('No solver present, unable to create shortcuts')
except AttributeError:
    print('No solver present, unable to create shortcuts')    

from .CBRead import readSBML3FBC, readSBML2FBA, readCOBRASBML
from .CBWrite import writeSBML3FBC, writeFVAtoCSV, writeModelToExcel97, writeModelToCOMBINEarchive, writeCOBRASBML, writeSBML3FBCV2


## "optional" modules extended dependencies
# debug
#import pyscescbm.CBGUI

try:
    from .CBGUI import *
except ImportError as ex:
    print(ex)
    print('GUI module not available')
    CBGUI = None
#import CBMultiEnv
try:
    from . import CBMultiEnv
except ImportError as ex:
    print(ex)
    print('MultiEnvironment module not available')
    CBMultiEnv = None
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

from . import CBModel, CBDataStruct, CBModelTools, CBRead, CBReadtxt, CBTools, CBVersion, CBWrite, CBXML, CBNetDB, CBPlot

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
    #print('*                     Somewhere In Time                               *')
    print('* Copyright(C) Brett G. Olivier 2010 - 2015                           *')
    print('* Dept. of Systems Bioinformatics                                     *')
    print('* Vrije Universiteit Amsterdam, Amsterdam, The Netherlands            *')
    #print('* PySCeS is distributed under the PySCeS (BSD style) licence, see     *')
    print('* CBMPy is distributed under the GNU GPL v 3.0 licence, see           *')
    print('* LICENCE (supplied with this release) for details                    *')
    ##  print('*                 ** Read about PySCeS **                             *')
    ##  print('* http://bioinformatics.oupjournals.org/cgi/content/short/21/4/560    *')
    print('***********************************************************************\n')


print('\nPySCeS-CBM is now CBMPy please update your scripts for future compatability.\nimport pyscescbm as cbm --> import cbmpy as cbm\n')