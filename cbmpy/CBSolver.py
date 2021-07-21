"""
CBMPy: CBSolver module
======================
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
Last edit: $Author: bgoli $ ($Id: CBSolver.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']
SOLVER_PREF = __CBCONFIG__['SOLVER_PREF']

__AVAILABLE_SOLVERS__ = []

__CBCONFIG__['SOLVER_ACTIVE'] = None
try:
    # for release 0.8.2+
    # from . import CBGLPK
    # from .CBGLPK import *

    # replaced by new module
    from . import CBGLPKOLD
    from .CBGLPKOLD import *

    __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK'
    __AVAILABLE_SOLVERS__.append('GLPK')
except ImportError:
    print('GLPK not available')
try:
    from . import CBCPLEX
    from .CBCPLEX import *

    if __CBCONFIG__['SOLVER_ACTIVE'] == None:
        __CBCONFIG__['SOLVER_ACTIVE'] = 'CPLEX'
    else:
        __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK+CPLEX'
    __AVAILABLE_SOLVERS__.append('CPLEX')
except Exception as ex:
    print('\n\n')
    print(ex)
    print('\n\n')
    print('CPLEX not available')

if __CBCONFIG__['SOLVER_ACTIVE'] == None:
    print(
        '\n*****\nWARNING: No linear solver present, please install IBM CPLEX with Python bindings or PyGLPK, please see http://cbmpy.sourceforge.net for Windows binary or http://tfinley.net/software/pyglpk for source.\n*****\n'
    )

if __CBCONFIG__['SOLVER_ACTIVE'] == 'GLPK+CPLEX':
    if __CBCONFIG__['SOLVER_PREF'] == None:
        __CBCONFIG__['SOLVER_ACTIVE'] = None
        while __CBCONFIG__['SOLVER_ACTIVE'] not in ['CPLEX', 'GLPK']:
            __CBCONFIG__['SOLVER_ACTIVE'] = raw_input(
                '\nSolver preference not set. Please select a solver [CPLEX or GLPK]: '
            )

    elif __CBCONFIG__['SOLVER_PREF'] == 'CPLEX':
        __CBCONFIG__['SOLVER_ACTIVE'] = 'CPLEX'
    elif __CBCONFIG__['SOLVER_PREF'] == 'GLPK':
        __CBCONFIG__['SOLVER_ACTIVE'] = 'GLPK'


__GLPK_METHODS__ = []
__CPLEX_METHODS__ = []
for k in tuple(globals()):
    if k[:5] == 'glpk_':
        __GLPK_METHODS__.append(k)
    elif k[:5] == 'cplx_':
        __CPLEX_METHODS__.append(k)

__COMMON_METHODS__ = [
    'analyzeModel',
    'FluxVariabilityAnalysis',
    'getOptimalSolution',
    'MinimizeSumOfAbsFluxes',
]

# print '__GLPK_METHODS__', __GLPK_METHODS__
# print '__CPLEX_METHODS__', __CPLEX_METHODS__
# print '__PUBLIC_METHODS__', __COMMON_METHODS__


def __setSolverInit__(slv):
    """
    Sets the active solver:

     - *slv* is either 'GLPK' or 'CPLEX'

    """
    if slv == 'GLPK':
        for k in tuple(globals()):
            if k[:5] == 'glpk_' and k[5:] in __COMMON_METHODS__:
                globals().update({k[5:]: globals()[k]})
    elif slv == 'CPLEX':
        for k in tuple(globals()):
            if k[:5] == 'cplx_' and k[5:] in __COMMON_METHODS__:
                globals().update({k[5:]: globals()[k]})


if __CBCONFIG__['SOLVER_ACTIVE'] == 'CPLEX':
    __setSolverInit__('CPLEX')
    print('\n*****\nUsing CPLEX\n*****\n')
elif __CBCONFIG__['SOLVER_ACTIVE'] == 'GLPK':
    from .CBGLPK import *

    __setSolverInit__('GLPK')
    print('\n*****\nUsing GLPK\n*****\n')
