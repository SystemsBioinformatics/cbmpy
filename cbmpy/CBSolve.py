"""
CBMPy: CBSolve module
=====================
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
Last edit: $Author: bgoli $ ($Id: CBModel.py 706 2020-03-23 21:31:49Z bgoli $)

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

import os
import time
import gc

# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio

import numpy
from . import CBWrite, CBTools
from .CBConfig import __CBCONFIG__ as __CBCONFIG__


from packaging import version as pkgver

HAVE_SYMPY = False
try:
    import sympy

    if pkgver.parse(sympy.__version__) >= pkgver.Version('0.7.5'):
        HAVE_SYMPY = True

    else:
        del sympy
        print(
            '\nWARNING: SymPy version 0.7.5 or newer is required for symbolic matrix support.'
        )
except ImportError:
    HAVE_SYMPY = False

HAVE_SCIPY = False
try:
    from scipy.sparse import csr

    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


class LPproperties:
    """LP property class"""

    __m = None
    __id = None
    __name = None

    @property
    def model(self):
        return self.__m

    @model.setter
    def model(self, model):
        self.__m = model

    @property
    def id(self):
        return self.__pid

    @id.setter
    def id(self, pid):
        self.__pid = pid

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, name):
        self.__name = name


class LPbase(LPproperties):
    """LP base class"""

    def __init__(self, pid):
        self.id = pid


class LPcplex(LPbase):
    """CPLEX interface"""

    def __init__(self, pid):
        super().__init__(pid)


class LPglpk(LPbase):
    """GLPK interface"""

    def __init__(self, pid):
        super().__init__(pid)


if __name__ == '__main__':
    lp = LPbase('base')
    lpc = LPcplex('cplex')
    lpg = LPglpk('glpk')
