"""
CBMPy: CBConfig module
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
Last edit: $Author: bgoli $ ($Id: CBConfig.py 711 2020-04-27 14:22:34Z bgoli $)

"""
# gets rid of "invalid variable name" info
# pylint: disable=C0103
# gets rid of "line to long" info
# pylint: disable=C0301
# use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

import platform

__VERSION_MAJOR__ = 0
__VERSION_MINOR__ = 8
__VERSION_MICRO__ = 4

__CBCONFIG__ = {
    'VERSION_MAJOR': __VERSION_MAJOR__,
    'VERSION_MINOR': __VERSION_MINOR__,
    'VERSION_MICRO': __VERSION_MICRO__,
    'VERSION_STATUS': '',
    'VERSION': '{}.{}.{}'.format(
        __VERSION_MAJOR__, __VERSION_MINOR__, __VERSION_MICRO__
    ),
    'DEBUG': False,
    'SOLVER_PREF': 'CPLEX',
    #'SOLVER_PREF': 'GLPK',
    'SOLVER_ACTIVE': None,
    'REVERSIBLE_SYMBOL': '<==>',
    'IRREVERSIBLE_SYMBOL': '-->',
    'HAVE_SBML2': False,
    'HAVE_SBML3': False,
    'CBMPY_DIR': None,
    'SYMPY_DENOM_LIMIT': 10 ** 32,
    'ENVIRONMENT': '{} {} ({})'.format(
        platform.system(), platform.release(), platform.architecture()[0]
    ),
    'MULTICORE_PYTHON_BIN_OVERRIDE': None,
    
}


def current_version():
    """
    Return the current CBMPy version as a string

    """
    return '{}.{}.{}'.format(__VERSION_MAJOR__, __VERSION_MINOR__, __VERSION_MICRO__)


def current_version_tuple():
    """
    Return the current CBMPy version as a tuple (x, y, z)

    """
    return (__VERSION_MAJOR__, __VERSION_MINOR__, __VERSION_MICRO__)
