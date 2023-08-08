"""
CBMPy: Constraint Based Modelling in Python (http://pysces.sourceforge.net/cbm)
============
Copyright (C) 2010-2024 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: setupegg.py 660 2024-09-24 14:57:04Z bgoli $)

"""

"""
A setup.py script to use setuptools, which gives egg goodness, etc.

Adapted from the original NumPy src (numpy.scipy.org).
"""
FRYING_EGGS = True
from setuptools import setup

execfile('setup.py')
