"""
CBMPy: __init__.py
==================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2022 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: __init__.py 710 2020-04-27 14:22:34Z bgoli $)

"""
#
# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

import os

# Setup some environment
__SILENT_START__ = False
from . import CBConfig


__CBCONFIG__ = CBConfig.__CBCONFIG__
__version__ = CBConfig.__CBCONFIG__['VERSION']
# Force the use of the GLPK solver by defining this envar

# Force the use of the CPLX solver by defining this envar
if 'CBMPY_USE_CPLX' in os.environ:
    CBConfig.__CBCONFIG__['SOLVER_PREF'] = 'CPLX'
if 'CBMPY_USE_GLPK' in os.environ:
    CBConfig.__CBCONFIG__['SOLVER_PREF'] = 'GLPK'

__CBCONFIG__['CBMPY_DIR'] = os.path.split(CBConfig.__file__)[0]

# This is just a hack for backwards compatibility with existing scripts
# that enables "from cbmpy import *" to work - bgoli

# lets add some shortcuts that will be core commands
__HAVE_SOLVER__ = False
try:
    from . import CBSolver

    analyzeModel = CBSolver.analyzeModel
    FluxVariabilityAnalysis = CBSolver.FluxVariabilityAnalysis

    # temporary hack
    if CBConfig.__CBCONFIG__['SOLVER_PREF'] == 'CPLX':
        MinimizeSumOfAbsFluxes = CBSolver.MinimizeSumOfAbsFluxes
        doFBAMinSum = CBSolver.MinimizeSumOfAbsFluxes
    else:
        print('doFBAMinSum not available with GLPK')

        def MinimizeSumOfAbsFluxes(*args, **kwargs):
            raise RuntimeError('doFBAMinSum not available with GLPK')

        doFBAMinSum = MinimizeSumOfAbsFluxes

    doFBA = CBSolver.analyzeModel
    doFVA = CBSolver.FluxVariabilityAnalysis
    # doFBAMinSum = CBSolver.MinimizeSumOfAbsFluxes

    __HAVE_SOLVER__ = True
except (ImportError, AttributeError) as e:
    print(e)
    print('No solver present, unable to create shortcuts')

from .CBWrite import writeFVAtoCSV, writeModelToExcel97
from .CBXML import _HAVE_SBML_

if _HAVE_SBML_:
    from .CBXML import sbml_fileFindVersion, sbml_fileValidate
    from .CBRead import loadModel, readSBML3FBC, readSBML2FBA, readCOBRASBML
    from .CBWrite import (
        saveModel,
        writeSBML3FBC,
        writeModelToExcel97,
        writeModelToCOMBINEarchive,
        writeCOBRASBML,
        writeSBML3FBCV2,
    )
else:
    print('\nWARNING: No SBML support, top-level SBML read/write functions disabled.\n')

try:
    from . import CBMultiCore
except ImportError as ex:
    print(ex)
    print('MultiCore module not available')
    CBMultiCore = None

INF = float('inf')
NINF = -float('inf')
NAN = float('nan')

from . import (
    CBModel,
    CBCommon,
    CBModelTools,
    CBRead,
    CBReadtxt,
    CBTools,
    CBWrite,
    CBXML,
    CBNetDB,
    CBPlot,
    CBMultiModel,
)
from .CBTools import deSerialize, deSerializeFromDisk
from . import PyscesSED

SED = PyscesSED.SEDCBMPY


def test():
    print(
        '\nCBMPy test suite is now separate from CMBPy source, see https://github.com/SystemsBioinformatics/cbmpy for details.\n'
    )


del os

if not __SILENT_START__:
    print('\n\n***********************************************************************')
    print(
        '* Welcome to CBMPy ({}) - PySCeS Constraint Based Modelling        *'.format(
            CBConfig.__CBCONFIG__['VERSION']
        )
    )
    print('*                http://cbmpy.sourceforge.net                         *')
    print('* Copyright(C) Brett G. Olivier 2014 - 2020                           *')
    print('* Systems Biology Lab, Vrije Universiteit Amsterdam                   *')
    print('* Amsterdam, The Netherlands                                          *')
    print('* CBMPy is developed as part of the BeBasic MetaToolKit Project       *')
    print('* Distributed under the GNU GPL v 3.0 licence, see                    *')
    print('* LICENCE (supplied with this release) for details                    *')
    print('***********************************************************************\n')
