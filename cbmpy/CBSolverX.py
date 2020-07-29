"""
CBMPy: CBSolver2 module
=======================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2018 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBSolver.py 305 2018-04-23 15:18:31Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']
SOLVER_PREF = __CBCONFIG__['SOLVER_PREF']
#__CBCONFIG__['SOLVER_ACTIVE'] = None

HAVE_CPLEX = False
HAVE_ECYGLPKI = False
HAVE_GUROBI = False
try:
    import cplex as cplex
    HAVE_CPLEX = True
except ImportError:
    pass
try:
    import ecyglpki as ecyglpki
    HAVE_ECYGLPKI = True
except ImportError as why:
    print(why)
# try:
    #import gurobi as _gurobi_
    #HAVE_GUROBI = True
# except ImportError:
    # pass
print('HAVE_ECYGLPKI:', HAVE_ECYGLPKI)
print('HAVE_CPLEX:', HAVE_CPLEX)
print('HAVE_GUROBI:', HAVE_GUROBI)


class LPBase(object):
    """
    The base of all LP's

    """
    obj = None
    id = None
    name = None
    objective = None
    objective_sense = None
    linear_constraints = None
    _constraint_senses_ = ('GE', 'LE', 'EQ')
    _objective_senses_ = ('MIN', 'MAX')

    def addObjective(self, arr):
        """
        Add objective as a list of (coefficient, variable) tuples

         - *arr* a list of one or more (coefficient, variable) tuples

        """
        self.objective = arr

    def setObjectiveSense(self, sense):
        """
        Set the objective sense.

         - *sense* is the string 'MAX' or 'MIN'

        """
        assert sense in self._objective_senses_, '\n{} is not a valid sense, must be one of {}.'.format(sense, self._objective_senses_)
        self.objective_sense = sense

    def addLinearConstraint(self, cid, arr, sense, rhs, overwrite=False):
        """
        Add linear constraint as a list of (coefficient, variable) tuples

         - *cid* the constraint id
         - *arr* a list of (coefficient, variable) tuples
         - *sense* GE, LE, EQ (>=, <=, ==)
         - *rhs* the floating point value of the rhs
         - *overwrite* [default=False] if a cid already exists overwrite it

        """
        assert len(arr) > 0, '\nAt least one coefficient, variable pair needs to be defined.'
        assert sense in self._constraint_senses_, '\n{} is not a valid sense, must be one of {}.'.format(sense, self._constraint_senses_)
        if self.linear_constraints is None:
            self.linear_constraints = {}
        if cid in self.linear_constraints and not overwrite:
            raise RuntimeWarning('\nConstraint id {} already exists.'.format(cid))
        else:
            self.linear_constraints[cid] = (arr, sense, rhs)


class MILPBase(object):
    """
    The base of all MILP's

    """


class QPbase(object):
    """
    The base of all QP's

    """


class LPCPLEX(LPBase):
        """
    CPLEX LP implimentation

    """


class LPECYGLPKI(LPBase):
        """
    ECYGLPKI LP implimentation

    """


c = LPCPLEX()
g = LPECYGLPKI()
