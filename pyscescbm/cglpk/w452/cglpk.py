# ctypes-glpk - A Python wrapper for GLPK using ctypes
# Copyright (C) 2008  Minh-Tri Pham

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# For further inquiries, please contact Minh-Tri Pham at pmtri80@gmail.com.
# ---------------------------------------------------------------------
#!/usr/bin/env python
"""ctypes-glpk - A ctypes-based Python wrapper for GLPK

ctypes-glpk is a Python module which encapsulates the functionality of the GNU Linear Programming Kit (GLPK). The GLPK allows one to specify linear programs (LPs) and mixed integer programs (MIPs), and to solve them with either simplex, interior-point, or branch-and-cut algorithms. The goal of ctypes-glpk is to give one Python access to all documented functionality of GLPK.

:Author: Minh-Tri Pham <pmtri80@gmail.com>
:Version: 0.2.3 (stable)
:Released: October 2008 (stable)

Availability
============

   To get the lastest version, see:
   http://code.google.com/p/ctypes-glpk/

Requirements
============

There are two requirements. First, GLPK must be installed on your platform. ctypes-glpk detects the existence of GLPK by calling 'glpsol -v'. It uses the corresponding shared library of GLPK to provide the functionality. Second, you need package 'ctypes' installed on your Python, which is a built-in starting from Python 2.5.

Installation
============

No installation is required. Just import this module in your code and run. You can also include this module in your <Python>\lib\site-packages folder, or include the path to the file in your PYTHONPATH variable.

Support
-------

- Platforms: Win32, Linux, and Mac OS
- GLPK: 4.9 to 4.33

How to use
==========

See the 'sample.c' file for an example of using GLPK in C, and the 'sample.py' file for the corresponding code in Python. You will find that they are almost identical.


Change Log
==========

ctypes-glpk-0.2.4 release
-------------------------

- Added wrapping functions for GLPK version 4.33

ctypes-glpk-0.2.3 release
-------------------------

- Fixed a number of bugs -- thanks Steve Jackson for having raised the issues

ctypes-glpk-0.2.2 release
-------------------------

- Fixed a bug in reading GLPK's version from glpsol -- thanks Steve Jackson for having raised the issue
- Fixed a bug in defining glp_erase_prob() -- thanks Steve Jackson for having raised the issue and having given a bug fix

ctypes-glpk-0.2.1 stable release
--------------------------------

- Added wrapping functions for GLPK version 4.32

ctypes-glpk-0.2.0 stable release
--------------------------------

- Added wrapping functions for GLPK from version 4.10 to version 4.31

ctypes-glpk-0.1.0 development release
-------------------------------------

- Added wrapping functions for GLPK version 4.9
- Created project ctypes-glpk

"""

#=============================================================================
# Importation
#=============================================================================
from ctypes import *
import os, sys

# Common declarations
c_int_p = POINTER(c_int)
c_float_p = POINTER(c_float)
c_double_p = POINTER(c_double)
size_t = c_uint


#=============================================================================
# Detect GLPK and load it
#=============================================================================
def _load_glpk():
    # Attempt to automatically detect which version of GLPK we are using.
    try:
        fi, fo = os.popen4('glpsol -v')
        fi.close()
        tokens = fo.read().split()

        # Version GLPK detected!!
        i = tokens.index('LP/MIP')
        version_string = tokens[i+3] if tokens[i+1] == 'Solver,' else tokens[i+2] # not safe, but works for now
        version = version_string.split('.')
        version = tuple([int(i) for i in version[:2]])

    except Exception as e:
        raise ImportError("Failed to run 'glpsol' to extract version number. GLPK may not be properly installed.")

    def _load_lib(dllname):
        try:
            return cdll.LoadLibrary(dllname)
        except:
            raise ImportError("Cannot import GLPK's shared library (" + dllname + "). Make sure its path is included in your system's PATH variable.")

    # Attempt to load the DLL
    if os.name == 'posix' and sys.platform.startswith('linux'):
        try:
            glpk_lib = _load_lib('libglpk.so')
        except ImportError:
            glpk_lib = _load_lib('libglpk.so.0')
    elif os.name == 'posix' and sys.platform.startswith('darwin'):
        glpk_lib = _load_lib('libglpk.dylib')
    elif os.name == 'nt':
        glpk_lib = _load_lib('glpk'+str(version[0])+str(version[1])+'.dll')
    else:
        raise ImportError('Platform '+str(os.name)+' is currently not supported.')

    return version, glpk_lib

_version, _glpk_lib = _load_glpk()

# return the current GLPK version
if _version >= (4, 9):
    def lpx_version():
        return str(_version[0])+'.'+str(_version[1])

#=============================================================================
# make function prototypes a bit easier to declare
#=============================================================================
def cfunc(name, result, *args):
    '''build and apply a ctypes prototype complete with parameter flags
    e.g.
cvMinMaxLoc = cfunc('cvMinMaxLoc', None,
                    ('image', POINTER(IplImage), 1),
                    ('min_val', POINTER(double), 2),
                    ('max_val', POINTER(double), 2),
                    ('min_loc', POINTER(CvPoint), 2),
                    ('max_loc', POINTER(CvPoint), 2),
                    ('mask', POINTER(IplImage), 1, None))
means locate cvMinMaxLoc in dll _glpk_lib, it returns nothing.
The first argument is an input image. The next 4 arguments are output, and the last argument is
input with an optional value. A typical call might look like:

min_val,max_val,min_loc,max_loc = cvMinMaxLoc(img)
    '''
    atypes = []
    aflags = []
    for arg in args:
        atypes.append(arg[1])
        aflags.append((arg[2], arg[0]) + arg[3:])
    return CFUNCTYPE(result, *atypes)((name, _glpk_lib), tuple(aflags))


#=============================================================================
# Identifiers
#=============================================================================

if _version >= (4, 9) and _version <= (4, 15):
    _lpx = 'glp_lpx_'
elif _version >= (4, 16):
    _lpx = '_glp_lpx_'
    _glp = 'glp_'


#=============================================================================
# Problem Object
#=============================================================================

if _version >= (4, 9) and _version <= (4, 15):
    def LPX_FIELDS():
        return [] # to be expanded in the future

    class LPX(Structure):
        _fields_ = LPX_FIELDS()

if _version >= (4, 16):
    def GLP_PROB_FIELDS():
        return [] # to be expanded in the future

    class glp_prob(Structure):
        _fields_ = GLP_PROB_FIELDS()

    LPX = glp_prob


#=============================================================================
# Karush-Kuhn-Tucker
#=============================================================================

if _version >= (4, 9):
    def LPXKKT_FIELDS():
        return [
            ('pe_ae_max', c_double), # largest absolute error
            ('pe_ae_row', c_int), # number of row with largest absolute error
            ('pe_re_max', c_double), # largest relative error
            ('pe_re_row', c_int), # number of row with largest relative error
            ('pe_quality', c_int), # quality of primal solution

            ('pb_ae_max', c_double), # largest absolute error
            ('pb_ae_ind', c_int), # number of variable with largest absolute error
            ('pb_re_max', c_double), # largest relative error
            ('pb_re_ind', c_int), # number of variable with largest relative error
            ('pb_quality', c_int), # quality of primal feasibility

            ('de_ae_max', c_double), # largest absolute error
            ('de_ae_col', c_int), # number of column with largest absolute error
            ('de_re_max', c_double), # largest relative error
            ('de_re_col', c_int), # number of column with largest relative error
            ('de_quality', c_int), # quality of dual solution

            ('db_ae_max', c_double), # largest absolute error
            ('db_ae_ind', c_int), # number of variable with largest absolute error
            ('db_re_max', c_double), # largest relative error
            ('db_re_ind', c_int), # number of variable with largest relative error
            ('db_quality', c_int), # quality of dual feasibility

            ('cs_ae_max', c_double), # largest absolute error
            ('cs_ae_ind', c_int), # number of variable with largest absolute error
            ('cs_re_max', c_double), # largest relative error
            ('cs_re_ind', c_int), # number of variable with largest relative error
            ('cs_quality', c_int), # quality of complementary slackness
        ]

    class LPXKKT(Structure):
        _fields_ = LPXKKT_FIELDS()


#=============================================================================
# Problem creating and modifying routines
#=============================================================================

if _version >= (4, 9):
    # create problem object
    lpx_create_prob = cfunc(_lpx+'create_prob', POINTER(LPX))

    # assign (change) problem name
    lpx_set_prob_name = cfunc(_lpx+'set_prob_name', None,
        ('lp', POINTER(LPX), 1),
        ('name', c_char_p, 1),
    )

    # assign (change) objective function name
    lpx_set_obj_name = cfunc(_lpx+'set_obj_name', None,
        ('lp', POINTER(LPX), 1),
        ('name', c_char_p, 1),
    )

    # set (change) optimization direction flag
    lpx_set_obj_dir = cfunc(_lpx+'set_obj_dir', None,
        ('lp', POINTER(LPX), 1),
        ('dir', c_int, 1),
    )
    LPX_MIN = 120 # minimization
    LPX_MAX = 121 # maximization

    # add new rows to problem object
    lpx_add_rows = cfunc(_lpx+'add_rows', c_int,
        ('lp', POINTER(LPX), 1),
        ('nrs', c_int, 1),
    )

    # add new columns to problem object
    lpx_add_cols = cfunc(_lpx+'add_cols', c_int,
        ('lp', POINTER(LPX), 1),
        ('ncs', c_int, 1),
    )

    # assign (change) row name
    lpx_set_row_name = cfunc(_lpx+'set_row_name', None,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
        ('name', c_char_p, 1),
    )

    # assign (change) column name
    lpx_set_col_name = cfunc(_lpx+'set_col_name', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('name', c_char_p, 1),
    )

    # set (change) row bounds
    lpx_set_row_bnds = cfunc(_lpx+'set_row_bnds', None,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
        ('type', c_int, 1),
        ('lb', c_double, 1),
        ('ub', c_double, 1),
    )
    LPX_FR = 110 #free variable
    LPX_LO = 111 # variable with lower bound
    LPX_UP = 112 # variable with upper bound
    LPX_DB = 113 # double-bounded variable
    LPX_FX = 114 # fixed variable

    # set (change) column bounds
    lpx_set_col_bnds = cfunc(_lpx+'set_col_bnds', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('type', c_int, 1),
        ('lb', c_double, 1),
        ('ub', c_double, 1),
    )

    # set (change) objective coefficient or constant term
    lpx_set_obj_coef = cfunc(_lpx+'set_obj_coef', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('coef', c_double, 1),
    )

    # set (replace) row of the constraint matrix
    lpx_set_mat_row = cfunc(_lpx+'set_mat_row', None,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # set (replace) column of the constraint matrix
    lpx_set_mat_col = cfunc(_lpx+'set_mat_col', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # load (replace) the whole constraint matrix
    lpx_load_matrix = cfunc(_lpx+'load_matrix', None,
        ('lp', POINTER(LPX), 1),
        ('ne', c_int, 1),
        ('ia', c_int_p, 1),
        ('ja', c_int_p, 1),
        ('ar', c_double_p, 1),
    )

    # delete rows from problem object
    lpx_del_rows = cfunc(_lpx+'del_rows', None,
        ('lp', POINTER(LPX), 1),
        ('nrs', c_int, 1),
        ('num', c_int_p, 1),
    )

    # delete columns from problem object
    lpx_del_cols = cfunc(_lpx+'del_cols', None,
        ('lp', POINTER(LPX), 1),
        ('ncs', c_int, 1),
        ('num', c_int_p, 1),
    )

# -- Begin of bug fix by Steve Jackson, Oct'08
if _version >= (4, 29):
    # Erase problem object content
    glp_erase_prob = cfunc(_glp+'erase_prob', None,
        ('lp', POINTER(glp_prob), 1),
    )
# -- End of bug fix by Steve Jackson, Oct'08

if _version >= (4, 9):
    # delete problem object
    lpx_delete_prob = cfunc(_lpx+'delete_prob', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 16):
    # create problem object
    glp_create_prob = cfunc(_glp+'create_prob', POINTER(glp_prob))

    # assign (change) problem name
    glp_set_prob_name = cfunc(_glp+'set_prob_name', None,
        ('lp', POINTER(glp_prob), 1),
        ('name', c_char_p, 1),
    )

    # assign (change) objective function name
    glp_set_obj_name = cfunc(_glp+'set_obj_name', None,
        ('lp', POINTER(glp_prob), 1),
        ('name', c_char_p, 1),
    )

    # set (change) optimization direction flag
    glp_set_obj_dir = cfunc(_glp+'set_obj_dir', None,
        ('lp', POINTER(glp_prob), 1),
        ('dir', c_int, 1),
    )
    GLP_MIN   = 1  # minimization
    GLP_MAX   = 2  # maximization

    # add new rows to problem object
    glp_add_rows = cfunc(_glp+'add_rows', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('nrs', c_int, 1),
    )

    # add new columns to problem object
    glp_add_cols = cfunc(_glp+'add_cols', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('ncs', c_int, 1),
    )

    # assign (change) row name
    glp_set_row_name = cfunc(_glp+'set_row_name', None,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('name', c_char_p, 1),
    )

    # assign (change) column name
    glp_set_col_name = cfunc(_glp+'set_col_name', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('name', c_char_p, 1),
    )

    # set (change) row bounds
    glp_set_row_bnds = cfunc(_glp+'set_row_bnds', None,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('type', c_int, 1),
        ('lb', c_double, 1),
        ('ub', c_double, 1),
    )
    GLP_FR = 1 #free variable
    GLP_LO = 2 # variable with lower bound
    GLP_UP = 3 # variable with upper bound
    GLP_DB = 4 # double-bounded variable
    GLP_FX = 5 # fixed variable

    # set (change) column bounds
    glp_set_col_bnds = cfunc(_glp+'set_col_bnds', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('type', c_int, 1),
        ('lb', c_double, 1),
        ('ub', c_double, 1),
    )

    # set (change) objective coefficient or constant term
    glp_set_obj_coef = cfunc(_glp+'set_obj_coef', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('coef', c_double, 1),
    )

    # set (replace) row of the constraint matrix
    glp_set_mat_row = cfunc(_glp+'set_mat_row', None,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # set (replace) column of the constraint matrix
    glp_set_mat_col = cfunc(_glp+'set_mat_col', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # load (replace) the whole constraint matrix
    glp_load_matrix = cfunc(_glp+'load_matrix', None,
        ('lp', POINTER(glp_prob), 1),
        ('ne', c_int, 1),
        ('ia', c_int_p, 1),
        ('ja', c_int_p, 1),
        ('ar', c_double_p, 1),
    )

    # delete rows from problem object
    glp_del_rows = cfunc(_glp+'del_rows', None,
        ('lp', POINTER(glp_prob), 1),
        ('nrs', c_int, 1),
        ('num', c_int_p, 1),
    )

    # delete columns from problem object
    glp_del_cols = cfunc(_glp+'del_cols', None,
        ('lp', POINTER(glp_prob), 1),
        ('ncs', c_int, 1),
        ('num', c_int_p, 1),
    )

if _version >= (4, 33):
    # copy problem object content
    glp_copy_prob = cfunc(_glp+'copy_prob', None,
        ('dest', POINTER(glp_prob), 1),
        ('prob', POINTER(glp_prob), 1),
        ('names', c_int, 1),
    )

if _version >= (4, 9):
    # delete problem object
    glp_delete_prob = cfunc(_glp+'delete_prob', None,
        ('lp', POINTER(glp_prob), 1),
    )


#=============================================================================
# Problem retrieving routines
#=============================================================================

if _version >= (4, 9):
    # retrieve problem name
    lpx_get_prob_name = cfunc(_lpx+'get_prob_name', c_char_p,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve objective function name
    lpx_get_obj_name = cfunc(_lpx+'get_obj_name', c_char_p,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve optimization direction flag
    lpx_get_obj_dir = cfunc(_lpx+'get_obj_dir', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve number of rows
    lpx_get_num_rows = cfunc(_lpx+'get_num_rows', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve number of columns
    lpx_get_num_cols = cfunc(_lpx+'get_num_cols', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve row name
    lpx_get_row_name = cfunc(_lpx+'get_row_name', c_char_p,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve column name
    lpx_get_col_name = cfunc(_lpx+'get_col_name', c_char_p,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve row type
    lpx_get_row_type = cfunc(_lpx+'get_row_type', c_int,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve row lower bound
    lpx_get_row_lb = cfunc(_lpx+'get_row_lb', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve row upper bound
    lpx_get_row_ub = cfunc(_lpx+'get_row_ub', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve column type
    lpx_get_col_type = cfunc(_lpx+'get_col_type', c_int,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve column lower bound
    lpx_get_col_lb = cfunc(_lpx+'get_col_lb', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve column upper bound
    lpx_get_col_ub = cfunc(_lpx+'get_col_ub', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve objective coefficient or constant term
    lpx_get_obj_coef = cfunc(_lpx+'get_obj_coef', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve number of constraint coefficients
    lpx_get_num_nz = cfunc(_lpx+'get_num_nz', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve row of the constraint matrix
    lpx_get_mat_row = cfunc(_lpx+'get_mat_row', c_int,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # retrieve column of the constraint matrix
    lpx_get_mat_col = cfunc(_lpx+'get_mat_col', c_int,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

if _version >= (4, 16):
    # retrieve problem name
    glp_get_prob_name = cfunc(_glp+'get_prob_name', c_char_p,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve objective function name
    glp_get_obj_name = cfunc(_glp+'get_obj_name', c_char_p,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve optimization direction flag
    glp_get_obj_dir = cfunc(_glp+'get_obj_dir', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve number of rows
    glp_get_num_rows = cfunc(_glp+'get_num_rows', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve number of columns
    glp_get_num_cols = cfunc(_glp+'get_num_cols', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve row name
    glp_get_row_name = cfunc(_glp+'get_row_name', c_char_p,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve column name
    glp_get_col_name = cfunc(_glp+'get_col_name', c_char_p,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve row type
    glp_get_row_type = cfunc(_glp+'get_row_type', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve row lower bound
    glp_get_row_lb = cfunc(_glp+'get_row_lb', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve row upper bound
    glp_get_row_ub = cfunc(_glp+'get_row_ub', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve column type
    glp_get_col_type = cfunc(_glp+'get_col_type', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve column lower bound
    glp_get_col_lb = cfunc(_glp+'get_col_lb', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve column upper bound
    glp_get_col_ub = cfunc(_glp+'get_col_ub', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve objective coefficient or constant term
    glp_get_obj_coef = cfunc(_glp+'get_obj_coef', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve number of constraint coefficients
    glp_get_num_nz = cfunc(_glp+'get_num_nz', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve row of the constraint matrix
    glp_get_mat_row = cfunc(_glp+'get_mat_row', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # retrieve column of the constraint matrix
    glp_get_mat_col = cfunc(_glp+'get_mat_col', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )


#=============================================================================
# Row and column searching routines
#=============================================================================

if _version >= (4, 15):
    # Create the name index
    lpx_create_index = cfunc(_lpx+'create_index', None,
        ('lp', POINTER(LPX), 1),
    )

    # Find row by its name
    lpx_find_row = cfunc(_lpx+'find_row', c_int,
        ('lp', POINTER(LPX), 1),
        ('name', c_char_p, 1),
    )

    # Find column by its name
    lpx_find_col = cfunc(_lpx+'find_col', c_int,
        ('lp', POINTER(LPX), 1),
        ('name', c_char_p, 1),
    )

    # Delete the name index
    lpx_delete_index = cfunc(_lpx+'delete_index', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 16):
    # Create the name index
    glp_create_index = cfunc(_glp+'create_index', None,
        ('lp', POINTER(glp_prob), 1),
    )

    # Find row by its name
    glp_find_row = cfunc(_glp+'find_row', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('name', c_char_p, 1),
    )

    # Find column by its name
    glp_find_col = cfunc(_glp+'find_col', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('name', c_char_p, 1),
    )

    # Delete the name index
    glp_delete_index = cfunc(_glp+'delete_index', None,
        ('lp', POINTER(glp_prob), 1),
    )


#=============================================================================
# Problem scaling routines
#=============================================================================

if _version >= (4, 18):
    # Set (change) row scale factor
    glp_set_rii = cfunc(_glp+'set_rii', None,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('rii', c_double, 1),
    )

    # Set (change) column scale factor
    glp_set_sjj = cfunc(_glp+'set_sjj', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('sjj', c_double, 1),
    )

    # Retrieve row scale factor
    glp_get_rii = cfunc(_glp+'get_rii', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # Retrieve column scale factor
    glp_get_sjj = cfunc(_glp+'get_sjj', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

if _version >= (4, 9):
    # scale problem data
    lpx_scale_prob = cfunc(_lpx+'scale_prob', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 31):
    # scaling options:
    GLP_SF_GM    = 0x01  # perform geometric mean scaling
    GLP_SF_EQ    = 0x10  # perform equilibration scaling
    GLP_SF_2N    = 0x20  # round scale factors to power of two
    GLP_SF_SKIP  = 0x40  # skip scaling if problem is well scaled
    GLP_SF_AUTO  = 0x80  # choose scaling options automatically

    # scale problem data
    glp_scale_prob = cfunc(_glp+'scale_prob', None,
        ('lp', POINTER(glp_prob), 1),
        ('flags', c_int, 1),
    )

if _version >= (4, 9):
    # unscale problem data
    lpx_unscale_prob = cfunc(_lpx+'unscale_prob', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 18):
    # unscale problem data
    glp_unscale_prob = cfunc(_glp+'unscale_prob', None,
        ('lp', POINTER(glp_prob), 1),
    )


#=============================================================================
# LP basis constructing routines
#=============================================================================

if _version >= (4, 9):
    # set (change) row status
    lpx_set_row_stat = cfunc(_lpx+'set_row_stat', None,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
        ('stat', c_int, 1),
    )
    LPX_BS = 140 # basic variable
    LPX_NL = 141 # non-basic variable on lower bound
    LPX_NU = 142 # non-basic variable on upper bound
    LPX_NF = 143 # non-basic free variable
    LPX_NS = 144 # non-basic fixed variable

    # set (change) column status
    lpx_set_col_stat = cfunc(_lpx+'set_col_stat', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('stat', c_int, 1),
    )

if _version >= (4, 16):
    # set (change) row status
    glp_set_row_stat = cfunc(_glp+'set_row_stat', None,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
        ('stat', c_int, 1),
    )
    GLP_BS = 1 # basic variable
    GLP_NL = 2 # non-basic variable on lower bound
    GLP_NU = 3 # non-basic variable on upper bound
    GLP_NF = 4 # non-basic free variable
    GLP_NS = 5 # non-basic fixed variable

    # set (change) column status
    glp_set_col_stat = cfunc(_glp+'set_col_stat', None,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('stat', c_int, 1),
    )

if _version >= (4, 9):
    # construct standard initial LP basis
    lpx_std_basis = cfunc(_lpx+'std_basis', None,
        ('lp', POINTER(LPX), 1),
    )

    # construct advanced initial LP basis
    lpx_adv_basis = cfunc(_lpx+'adv_basis', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 31):
    # construct standard initial LP basis
    glp_std_basis = cfunc(_glp+'std_basis', None,
        ('lp', POINTER(glp_prob), 1),
    )

    # construct advanced initial LP basis
    glp_adv_basis = cfunc(_glp+'adv_basis', None,
        ('lp', POINTER(glp_prob), 1),
        ('flags', c_int, 1),
    )

if _version >= (4, 10):
    # construct advanced initial LP basis
    lpx_cpx_basis = cfunc(_lpx+'cpx_basis', None,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 31):
    # Construct Bixby's initial LP basis
    glp_cpx_basis = cfunc(_glp+'cpx_basis', None,
        ('lp', POINTER(glp_prob), 1),
    )


#=============================================================================
# Simplex method routine
#=============================================================================

if _version >= (4, 9):
    # solve LP problem using the simplex method
    lpx_simplex = cfunc(_lpx+'simplex', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # exit codes returned by solver routines:
    LPX_E_OK = 200 # success
    LPX_E_EMPTY = 201 # empty problem
    LPX_E_BADB = 202 # invalid initial basis
    LPX_E_INFEAS = 203 # infeasible initial solution
    LPX_E_FAULT = 204 # unable to start the search
    LPX_E_OBJLL = 205 # objective lower limit reached
    LPX_E_OBJUL = 206 # objective upper limit reached
    LPX_E_ITLIM = 207 # iterations limit exhausted
    LPX_E_TMLIM = 208 # time limit exhausted
    LPX_E_NOFEAS = 209 # no feasible solution
    LPX_E_INSTAB = 210 # numerical instability
    LPX_E_SING = 211 # problems with basis matrix
    LPX_E_NOCONV = 212 # no convergence (interior)
    LPX_E_NOPFS = 213 # no primal feas. sol. (LP presolver)
    LPX_E_NODFS = 214 # no dual feas. sol. (LP presolver)

if _version >= (4, 32):
    LPX_E_MIPGAP    = 215   # relative mip gap tolerance reached

if _version >= (4, 18):
    def GLP_SMCP_FIELDS():
        return [ # simplex method control parameters
            ('msg_lev', c_int), # message level
            ('meth', c_int), # simplex method option
            ('pricing', c_int), # pricing technique
            ('r_test', c_int), # ratio test technique
            ('tol_bnd', c_double), # spx.tol_bnd
            ('tol_dj', c_double), # spx.tol_dj
            ('tol_piv', c_double), # spx.tol_piv
            ('obj_ll', c_double), # spx.obj_ll
            ('obj_ul', c_double), # spx.obj_ul
            ('it_lim', c_int), # spx.it_lim
            ('tm_lim', c_int), # spx.tm_lim (miliseconds)
            ('out_frq', c_int), # spx.out_frq
            ('out_dly', c_int), # spx.out_dly (miliseconds)
            ('presolve', c_int), # enable/disable using LP presolver
            ('foo_bar', c_double*36), # (reserved)
        ]

    class glp_smcp(Structure):
        _fields_ = GLP_SMCP_FIELDS()

    GLP_MSG_OFF     = 0  # no output
    GLP_MSG_ERR     = 1  # warning and error messages only
    GLP_MSG_ON      = 2  # normal output
    GLP_MSG_ALL     = 3  # full output
if _version >= (4, 20):
    GLP_MSG_DBG     = 4  # debug output

if _version >= (4, 18):
    GLP_PRIMAL      = 1  # use primal simplex
    GLP_DUALP       = 2  # use dual; if it fails, use primal

if _version >= (4, 31):
    GLP_DUAL        = 3  # use dual simplex

if _version >= (4, 18):
    GLP_PT_STD   = 0x11  # standard (textbook)
    GLP_PT_PSE   = 0x22  # projected steepest edge

    GLP_RT_STD   = 0x11  # standard (textbook)
    GLP_RT_HAR   = 0x22  # Harris' ratio test

    # enable/disable flag:
    GLP_ON          = 1  # enable something
    GLP_OFF         = 0  # disable something

    # solve LP problem using the simplex method
    glp_simplex = cfunc(_glp+'simplex', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('parm', POINTER(glp_smcp), 1),
    )

    # return codes:
    GLP_EBADB    = 0x01  # invalid basis
    GLP_ESING    = 0x02  # singular matrix
    GLP_ECOND    = 0x03  # ill-conditioned matrix
    GLP_EBOUND   = 0x04  # invalid bounds
    GLP_EFAIL    = 0x05  # solver failed
    GLP_EOBJLL   = 0x06  # objective lower limit reached
    GLP_EOBJUL   = 0x07  # objective upper limit reached
    GLP_EITLIM   = 0x08  # iteration limit exceeded
    GLP_ETMLIM   = 0x09  # time limit exceeded
    GLP_ENOPFS   = 0x0A  # no primal feasible solution
    GLP_ENODFS   = 0x0B  # no dual feasible solution

if _version >= (4, 20):
    GLP_EROOT    = 0x0C  # root LP optimum not provided
    GLP_ESTOP    = 0x0D  # serach terminated by application

if _version >= (4, 32):
    GLP_EMIPGAP  = 0x0E  # relative mip gap tolerance reached

if _version >= (4, 33):
    GLP_ENOFEAS     = 0x0F  # no primal/dual feasible solution
    GLP_ENOCVG      = 0x10  # no convergence
    GLP_EINSTAB     = 0x11  # numerical instability


if _version >= (4, 18):
    # Initialize simplex method control parameters
    glp_init_smcp = cfunc(_glp+'init_smcp', None,
        ('parm', POINTER(glp_smcp), 1),
    )
    # Return type on GLPK doc is int, while the actual return type is None -- found and fixed by Steve Jackson

if _version >= (4, 13):
    # solve LP problem using the primal two-phase simplex method based on exact (rational) arithmetic
    lpx_exact = cfunc(_lpx+'exact', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 33):
    # solve LP problem in exact arithmetic
    glp_exact = cfunc(_glp+'exact', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('parm', POINTER(glp_smcp), 1),
    )


#=============================================================================
# Basic solution retrieving routines
#=============================================================================

if _version >= (4, 9):
    # retrieve generic status of basic solution
    lpx_get_status = cfunc(_lpx+'get_status', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # status codes reported by the routine lpx_get_status:
    LPX_OPT = 180 # optimal
    LPX_FEAS = 181 # feasible
    LPX_INFEAS = 182 # infeasible
    LPX_NOFEAS = 183 # no feasible
    LPX_UNBND = 184 # unbounded
    LPX_UNDEF = 185 # undefined

    # retrieve primal status of basic solution
    lpx_get_prim_stat = cfunc(_lpx+'get_prim_stat', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # status of primal basic solution:
    LPX_P_UNDEF = 132 # primal solution is undefined
    LPX_P_FEAS = 133 # solution is primal feasible
    LPX_P_INFEAS = 134 # solution is primal infeasible
    LPX_P_NOFEAS = 135 # no primal feasible solution exists

    # retrieve dual status of basic solution
    lpx_get_dual_stat = cfunc(_lpx+'get_dual_stat', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # status of dual basic solution:
    LPX_D_UNDEF = 136 # dual solution is undefined
    LPX_D_FEAS = 137 # solution is dual feasible
    LPX_D_INFEAS = 138 # solution is dual infeasible
    LPX_D_NOFEAS = 139 # no dual feasible solution exists

    # retrieve objective value
    lpx_get_obj_val = cfunc(_lpx+'get_obj_val', c_double,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve row status
    lpx_get_row_stat = cfunc(_lpx+'get_row_stat', c_int,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve row primal value
    lpx_get_row_prim = cfunc(_lpx+'get_row_prim', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve row dual value
    lpx_get_row_dual = cfunc(_lpx+'get_row_dual', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve column status
    lpx_get_col_stat = cfunc(_lpx+'get_col_stat', c_int,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve col primal value
    lpx_get_col_prim = cfunc(_lpx+'get_col_prim', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve column dual value
    lpx_get_col_dual = cfunc(_lpx+'get_col_dual', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

if _version >= (4, 16):
    # retrieve generic status of basic solution
    glp_get_status = cfunc(_glp+'get_status', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # status codes reported by the routine glp_get_status:
    GLP_UNDEF       = 1  # solution is undefined
    GLP_FEAS        = 2  # solution is feasible
    GLP_INFEAS      = 3  # solution is infeasible
    GLP_NOFEAS      = 4  # no feasible solution exists
    GLP_OPT         = 5  # solution is optimal
    GLP_UNBND       = 6  # solution is unbounded

    # retrieve primal status of basic solution
    glp_get_prim_stat = cfunc(_glp+'get_prim_stat', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve dual status of basic solution
    glp_get_dual_stat = cfunc(_glp+'get_dual_stat', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve objective value
    glp_get_obj_val = cfunc(_glp+'get_obj_val', c_double,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve row status
    glp_get_row_stat = cfunc(_glp+'get_row_stat', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve row primal value
    glp_get_row_prim = cfunc(_glp+'get_row_prim', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve row dual value
    glp_get_row_dual = cfunc(_glp+'get_row_dual', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve column status
    glp_get_col_stat = cfunc(_glp+'get_col_stat', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve col primal value
    glp_get_col_prim = cfunc(_glp+'get_col_prim', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve column dual value
    glp_get_col_dual = cfunc(_glp+'get_col_dual', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

if _version >= (4, 9):
    # retrieve non-basic variable which causes unboundness
    lpx_get_ray_info = cfunc(_lpx+'get_ray_info', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 33):
    # determine variable causing unboundness
    glp_get_unbnd_ray = cfunc(_glp+'get_unbnd_ray', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

if _version >= (4, 9):
    # check Karush-Kuhn-Tucker conditions
    lpx_check_kkt = cfunc(_lpx+'check_kkt', None,
        ('lp', POINTER(LPX), 1),
        ('scaled', c_int, 1),
        ('kkt', POINTER(LPXKKT), 1),
    )


#=============================================================================
# Interior-point method routines
#=============================================================================

if _version >= (4, 9):
    # solve LP problem using the primal-dual interior point method
    lpx_interior = cfunc(_lpx+'interior', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 33):
    # solve LP problem with the interior point method
    glp_interior = cfunc(_glp+'interior', c_int,
        ('lp', POINTER(LPX), 1),
        ('parm', c_void_p, 1),
    )

if _version >= (4, 9):
    # retrieve status of interior-point solution
    lpx_ipt_status = cfunc(_lpx+'ipt_status', c_int,
        ('lp', POINTER(LPX), 1),
    )
    LPX_T_UNDEF = 150 # interior solution is undefined
    LPX_T_OPT = 151 # interior solution is optimal

if _version >= (4, 18):
    # retrieve status of interior-point solution
    glp_ipt_status = cfunc(_glp+'ipt_status', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

if _version >= (4, 9):
    # retrieve objective value
    lpx_ipt_obj_val = cfunc(_lpx+'ipt_obj_val', c_double,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve row primal value
    lpx_ipt_row_prim = cfunc(_lpx+'ipt_row_prim', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve row dual value
    lpx_ipt_row_dual = cfunc(_lpx+'ipt_row_dual', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve column primal value
    lpx_ipt_col_prim = cfunc(_lpx+'ipt_col_prim', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve column dual value
    lpx_ipt_col_dual = cfunc(_lpx+'ipt_col_dual', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

if _version >= (4, 16):
    # retrieve objective value
    glp_ipt_obj_val = cfunc(_glp+'ipt_obj_val', c_double,
        ('lp', POINTER(glp_prob), 1),
    )

    # retrieve row primal value
    glp_ipt_row_prim = cfunc(_glp+'ipt_row_prim', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve row dual value
    glp_ipt_row_dual = cfunc(_glp+'ipt_row_dual', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve column primal value
    glp_ipt_col_prim = cfunc(_glp+'ipt_col_prim', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve column dual value
    glp_ipt_col_dual = cfunc(_glp+'ipt_col_dual', c_double,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )


#=============================================================================
# Mixed Integer Programming routines
#=============================================================================

if _version >= (4, 9):
    # set (change) problem class
    lpx_set_class = cfunc(_lpx+'set_class', None,
        ('lp', POINTER(LPX), 1),
        ('klass', c_int, 1),
    )
    LPX_LP = 100 # linear programming (LP)
    LPX_MIP = 101 # mixed integer programming (MIP)

    # retrieve problem class
    lpx_get_class = cfunc(_lpx+'get_class', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 9):
    # set (change) column kind
    lpx_set_col_kind = cfunc(_lpx+'set_col_kind', None,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
        ('kind', c_int, 1),
    )
    LPX_CV = 160 # continuous variable
    LPX_IV = 161 # integer variable

    # retrieve column kind
    lpx_get_col_kind = cfunc(_lpx+'get_col_kind', c_int,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

    # retrieve number of integer columns
    lpx_get_num_int = cfunc(_lpx+'get_num_int', c_int,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve number of binary columns
    lpx_get_num_bin = cfunc(_lpx+'get_num_bin', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 16):
    # set (change) column kind
    glp_set_col_kind = cfunc(_glp+'set_col_kind', None,
        ('mip', POINTER(glp_prob), 1),
        ('j', c_int, 1),
        ('kind', c_int, 1),
    )
    GLP_CV   = 1 # continuous variable
    GLP_IV   = 2 # integer variable
    GLP_BV   = 3  # binary variable

    # retrieve column kind
    glp_get_col_kind = cfunc(_glp+'get_col_kind', c_int,
        ('mip', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # retrieve number of integer columns
    glp_get_num_int = cfunc(_glp+'get_num_int', c_int,
        ('mip', POINTER(glp_prob), 1),
    )

    # retrieve number of binary columns
    glp_get_num_bin = cfunc(_glp+'get_num_bin', c_int,
        ('mip', POINTER(glp_prob), 1),
    )

if _version >= (4, 9):
    # solve MIP problem using the branch-and-bound method
    lpx_integer = cfunc(_lpx+'integer', c_int,
        ('mip', POINTER(LPX), 1),
    )

if _version >= (4, 20):
    # branch-and-bound tree
    class glp_tree(Structure):
        _fields_ = [('_tree', c_double),]

    glp_iocp_cback_func = CFUNCTYPE(None, glp_tree, c_void_p)

    def GLP_IOCP_FIELDS():
        return [ # integer optimizer control parameters
            ('msg_lev', c_int), # message level
            ('br_tech', c_int), # branching technique
            ('bt_tech', c_int), # backtracking technique
            ('tol_int', c_double), # mip.tol_int
            ('tol_obj', c_double), # mip.tol_obj
            ('tm_lim', c_int), # mip.tm_lim (miliseconds)
            ('out_frq', c_int), # mip.out_frq (miliseconds)
            ('out_dly', c_int), # mip.out_dly (miliseconds)
            ('cb_func', glp_iocp_cback_func), # mip.cb_func
            ('cb_info', c_void_p), # mip.cb_info
            ('cb_size', c_int), # mip.cb_size
            ('pp_tech', c_int), # preprocessing technique
            ('mip_gap', c_double), # relative MIP gap tolerance
            ('mir_cuts', c_int), # MIR cuts (GLP_ON/GLP_OFF)
            ('gmi_cuts', c_int), # Gomory's cuts (GLP_ON/GLP_OFF)
            ('cov_cuts', c_int), # cover cuts     (GLP_ON/GLP_OFF)
            ('clq_cuts', c_int), # clique cuts    (GLP_ON/GLP_OFF)
            ('presolve', c_int), # enable/disable using MIP presolver
            ('binarize', c_int), # try to binarize integer variables
            ('foo_bar', c_double*30), # (reserved)
            ('fn_sol', c_char_p), # file name to write solution found
        ]

    class glp_iocp(Structure):
        _fields_ = GLP_IOCP_FIELDS()

    GLP_BR_FFV      = 1  # first fractional variable
    GLP_BR_LFV      = 2  # last fractional variable
    GLP_BR_MFV      = 3  # most fractional variable
    GLP_BR_DTH      = 4  # heuristic by Driebeck and Tomlin

    GLP_BT_DFS      = 1  # depth first search
    GLP_BT_BFS      = 2  # breadth first search
    GLP_BT_BLB      = 3  # best local bound
    GLP_BT_BPH      = 4  # best projection heuristic

if _version >= (4, 21):
    GLP_PP_NONE     = 0  # disable preprocessing
    GLP_PP_ROOT     = 1  # preprocessing only on root level
    GLP_PP_ALL      = 2  # preprocessing on all levels

if _version >= (4, 20):
    # Solve MIP problem with the branch-and-cut method
    glp_intopt = cfunc(_glp+'intopt', c_int,
        ('mip', POINTER(glp_prob), 1),
        ('parm', POINTER(glp_iocp), 1),
    )

    # Initialize integer optimizer control parameters
    glp_init_iocp = cfunc(_glp+'init_iocp', None,
        ('parm', POINTER(glp_iocp), 1),
    )

if _version >= (4, 9):
    # solve MIP problem using the advanced B&B solver
    lpx_intopt = cfunc(_lpx+'intopt', c_int,
        ('mip', POINTER(LPX), 1),
    )

    # retrieve status of MIP solution
    lpx_mip_status = cfunc(_lpx+'mip_status', c_int,
        ('mip', POINTER(LPX), 1),
    )
    LPX_I_UNDEF = 170 # integer solution is undefined
    LPX_I_OPT = 171 # integer solution is optimal
    LPX_I_FEAS = 172 # integer solution is feasible
    LPX_I_NOFEAS = 173 # no integer solution exists

if _version >= (4, 18):
    # retrieve status of MIP solution
    glp_mip_status = cfunc(_glp+'mip_status', c_int,
        ('mip', POINTER(glp_prob), 1),
    )

if _version >= (4, 9):
    # retrieve objective value
    lpx_mip_obj_val = cfunc(_lpx+'mip_obj_val', c_double,
        ('lp', POINTER(LPX), 1),
    )

    # retrieve row value
    lpx_mip_row_val = cfunc(_lpx+'mip_row_val', c_double,
        ('lp', POINTER(LPX), 1),
        ('i', c_int, 1),
    )

    # retrieve column value
    lpx_mip_col_val = cfunc(_lpx+'mip_col_val', c_double,
        ('lp', POINTER(LPX), 1),
        ('j', c_int, 1),
    )

if _version >= (4, 16):
    # retrieve objective value
    glp_mip_obj_val = cfunc(_glp+'mip_obj_val', c_double,
        ('mip', POINTER(glp_prob), 1),
    )

    # retrieve row value
    glp_mip_row_val = cfunc(_glp+'mip_row_val', c_double,
        ('mip', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # retrieve column value
    glp_mip_col_val = cfunc(_glp+'mip_col_val', c_double,
        ('mip', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )


#=============================================================================
# Utility routines
#=============================================================================

if _version >= (4, 9):
    # read problem data in fixed MPS format
    lpx_read_mps = cfunc(_lpx+'read_mps', POINTER(LPX),
        ('fname', c_char_p, 1),
    )

    # write problem data in fixed MPS format
    lpx_write_mps = cfunc(_lpx+'write_mps', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # read LP basis in fixed MPS format
    lpx_read_bas = cfunc(_lpx+'read_bas', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # write LP basis in fixed MPS format
    lpx_write_bas = cfunc(_lpx+'write_bas', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # read problem data in free MPS format
    lpx_read_freemps = cfunc(_lpx+'read_freemps', POINTER(LPX),
        ('fname', c_char_p, 1),
    )

    # write problem data in free MPS format
    lpx_write_freemps = cfunc(_lpx+'write_freemps', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # read problem data in CLPEX LP format
    lpx_read_cpxlp = cfunc(_lpx+'read_cpxlp', POINTER(LPX),
        ('fname', c_char_p, 1),
    )

    # write problem data in CLPEX LP format
    lpx_write_cpxlp = cfunc(_lpx+'write_cpxlp', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

if _version >= (4, 29):
    # MPS file format:
    GLP_MPS_DECK    = 1  # fixed (ancient)
    GLP_MPS_FILE    = 2  # free (modern)

    # read problem data in MPS format
    glp_read_mps = cfunc(_glp+'read_mps', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fmt', c_int, 1),
        ('parm', c_void_p, 1),
        ('fname', c_char_p, 1),
    )

    # write problem data in MPS format
    glp_write_mps = cfunc(_glp+'write_mps', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fmt', c_int, 1),
        ('parm', c_void_p, 1),
        ('fname', c_char_p, 1),
    )

    # read problem data in CPLEX LP format
    glp_read_lp = cfunc(_glp+'read_lp', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('parm', c_void_p, 1),
        ('fname', c_char_p, 1),
    )

    # write problem data in CPLEX LP format
    glp_write_lp = cfunc(_glp+'write_lp', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('parm', c_void_p, 1),
        ('fname', c_char_p, 1),
    )

if _version >= (4, 9):
    # read model written in GNU MathProg modeling language
    lpx_read_model = cfunc(_lpx+'read_model', POINTER(LPX),
        ('model', c_char_p, 1),
        ('data', c_char_p, 1),
        ('output', c_char_p, 1),
    )

    # write problem data in plain text format
    lpx_print_prob = cfunc(_lpx+'print_prob', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # write basic solution in printable format
    lpx_print_sol = cfunc(_lpx+'print_sol', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # write bounds sensitivity information
    lpx_print_sens_bnds = cfunc(_lpx+'print_sens_bnds', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # write interior point solution in printable format
    lpx_print_ips = cfunc(_lpx+'print_ips', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

    # write MIP solution in printable format
    lpx_print_mip = cfunc(_lpx+'print_mip', c_int,
        ('lp', POINTER(LPX), 1),
        ('fname', c_char_p, 1),
    )

if _version >= (4, 23):
    # Read basic solution from text file
    glp_read_sol = cfunc(_glp+'read_sol', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )

    # Write basic solution to text file
    glp_write_sol = cfunc(_glp+'write_sol', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )

    # Read interior-point solution from text file
    glp_read_ipt = cfunc(_glp+'read_ipt', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )

    # Write interior-point solution to text file
    glp_write_ipt = cfunc(_glp+'write_ipt', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )

    # Read MIP solution from text file
    glp_read_mip = cfunc(_glp+'read_mip', c_int,
        ('mip', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )

    # Write MIP solution to text file
    glp_write_mip = cfunc(_glp+'write_mip', c_int,
        ('mip', POINTER(glp_prob), 1),
        ('fname', c_char_p, 1),
    )


#=============================================================================
# Subroutines for processing MathProg models
#=============================================================================

if _version >= (4, 33):
    # MathProg translator workspace
    class glp_tran(Structure):
        _fields_ = [('_opaque_tran', c_double),]

    # MathProg solution indicator:
    GLP_MPL_SOL        = 1  # basic solution
    GLP_MPL_IPT        = 2  # interior-point solution
    GLP_MPL_MIP        = 3  # mixed integer solution

    # allocate the MathProg translator workspace
    glp_mpl_alloc_wksp = cfunc(_glp+'mpl_alloc_wksp', POINTER(glp_tran),
    )

    # read and translate model section
    glp_mpl_read_model = cfunc(_glp+'mpl_read_model', c_int,
        ('tran', POINTER(glp_tran), 1),
        ('fname', c_char_p, 1),
        ('skip', c_int, 1),
    )

    # read and translate data section
    glp_mpl_read_data = cfunc(_glp+'mpl_read_data', c_int,
        ('tran', POINTER(glp_tran), 1),
        ('fname', c_char_p, 1),
    )

    # generate the model
    glp_mpl_generate = cfunc(_glp+'mpl_generate', c_int,
        ('tran', POINTER(glp_tran), 1),
        ('fname', c_char_p, 1),
    )

    # build LP/MIP problem instance from the model
    glp_mpl_build_prob = cfunc(_glp+'mpl_build_prob', None,
        ('tran', POINTER(glp_tran), 1),
        ('prob', POINTER(glp_prob), 1),
    )

    # postsolve the model
    glp_mpl_postsolve = cfunc(_glp+'mpl_postsolve', c_int,
        ('tran', POINTER(glp_tran), 1),
        ('prob', POINTER(glp_prob), 1),
        ('sol', c_int, 1),
    )

    # free the MathProg translator workspace
    glp_mpl_free_wksp = cfunc(_glp+'mpl_free_wksp', None,
        ('tran', POINTER(glp_tran), 1),
    )


#=============================================================================
# Control parameters and statistics routines
#=============================================================================

if _version >= (4, 9):
    # reset control parameters to default values
    lpx_reset_parms = cfunc(_lpx+'reset_parms', None,
        ('lp', POINTER(LPX), 1),
    )

    # set (change) integer control parameter
    lpx_set_int_parm = cfunc(_lpx+'set_int_parm', None,
        ('lp', POINTER(LPX), 1),
        ('parm', c_int, 1),
        ('val', c_int, 1),
    )

    # query integer control parameter
    lpx_get_int_parm = cfunc(_lpx+'get_int_parm', c_int,
        ('lp', POINTER(LPX), 1),
        ('parm', c_int, 1),
    )

    # set (change) real control parameter
    lpx_set_real_parm = cfunc(_lpx+'set_real_parm', None,
        ('lp', POINTER(LPX), 1),
        ('parm', c_int, 1),
        ('val', c_double, 1),
    )

    # query real control parameter
    lpx_get_real_parm = cfunc(_lpx+'get_real_parm', c_double,
        ('lp', POINTER(LPX), 1),
        ('parm', c_int, 1),
    )

    # control parameter identifiers:
    LPX_K_MSGLEV    = 300   # lp->msg_lev
    LPX_K_SCALE     = 301   # lp->scale
    LPX_K_DUAL      = 302   # lp->dual
    LPX_K_PRICE     = 303   # lp->price
    LPX_K_RELAX     = 304   # lp->relax
    LPX_K_TOLBND    = 305   # lp->tol_bnd
    LPX_K_TOLDJ     = 306   # lp->tol_dj
    LPX_K_TOLPIV    = 307   # lp->tol_piv
    LPX_K_ROUND     = 308   # lp->round
    LPX_K_OBJLL     = 309   # lp->obj_ll
    LPX_K_OBJUL     = 310   # lp->obj_ul
    LPX_K_ITLIM     = 311   # lp->it_lim
    LPX_K_ITCNT     = 312   # lp->it_cnt
    LPX_K_TMLIM     = 313   # lp->tm_lim
    LPX_K_OUTFRQ    = 314   # lp->out_frq
    LPX_K_OUTDLY    = 315   # lp->out_dly
    LPX_K_BRANCH    = 316   # lp->branch
    LPX_K_BTRACK    = 317   # lp->btrack
    LPX_K_TOLINT    = 318   # lp->tol_int
    LPX_K_TOLOBJ    = 319   # lp->tol_obj
    LPX_K_MPSINFO   = 320   # lp->mps_info
    LPX_K_MPSOBJ    = 321   # lp->mps_obj
    LPX_K_MPSORIG   = 322   # lp->mps_orig
    LPX_K_MPSWIDE   = 323   # lp->mps_wide
    LPX_K_MPSFREE   = 324   # lp->mps_free
    LPX_K_MPSSKIP   = 325   # lp->mps_skip
    LPX_K_LPTORIG   = 326   # lp->lpt_orig
    LPX_K_PRESOL    = 327   # lp->presol
    LPX_K_BINARIZE  = 328   # lp->binarize
    LPX_K_USECUTS   = 329   # lp->use_cuts

if _version >= (4, 17):
    LPX_K_BFTYPE    = 330   # lp->bf_type

if _version >= (4, 24):
    LPX_K_MIPGAP    = 331   # lp->mip_gap

if _version >= (4, 10):
    # control parameter identifiers:
    LPX_C_COVER     = 0x01  # mixed cover cuts
    LPX_C_CLIQUE    = 0x02  # clique cuts
    LPX_C_GOMORY    = 0x04  # Gomory's mixed integer cuts
    LPX_C_ALL       = 0xFF  # all cuts

if _version >= (4, 23):
    LPX_C_MIR       = 0x08  # mixed integer rounding cuts

#=============================================================================
# LP basis and simplex tableau routines
#=============================================================================

if _version >= (4, 18):
    # Check if the basis factorization exists
    glp_bf_exists = cfunc(_glp+'bf_exists', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # Compute the basis factorization
    glp_factorize = cfunc(_glp+'factorize', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    # Check if the basis factorization has been updated
    glp_bf_updated = cfunc(_glp+'bf_updated', c_int,
        ('lp', POINTER(glp_prob), 1),
    )

    def GLP_BFCP_FIELDS():
        return [ # basis factorization control parameters
            ('msg_lev', c_int), # message level
            ('type', c_int), # factorization type
            ('lu_size', c_int), # luf.sv_size
            ('piv_tol', c_double), # luf.piv_tol
            ('piv_lim', c_int), # luf.piv_lim
            ('suhl', c_int), # luf.suhl
            ('eps_tol', c_double), # luf.eps_tol
            ('max_gro', c_double), # luf.max_gro
            ('nfs_max', c_int), # fhv.hh_max
            ('upd_tol', c_double), # fhv.upd_tol
            ('nrs_max', c_int), # lpf.v_size
            ('rs_size', c_int), # lpf.rs_size
            ('foo_bar', c_double*38), # (reserved)
        ]

    GLP_BF_FT       = 1  # LUF + Forrest-Tomlin
    GLP_BF_BG       = 2  # LUF + Schur compl. + Bartels-Golub
    GLP_BF_GR       = 3  # LUF + Schur compl. + Givens rotation

    class glp_bfcp(Structure):
        _fields_ = GLP_BFCP_FIELDS()

    # Retrieve basis factorization control parameters
    glp_get_bfcp = cfunc(_glp+'get_bfcp', None,
        ('lp', POINTER(glp_prob), 1),
        ('parm', POINTER(glp_bfcp), 1),
    )

    # Change basis factorization control parameters
    glp_set_bfcp = cfunc(_glp+'set_bfcp', None,
        ('lp', POINTER(glp_prob), 1),
        ('parm', POINTER(glp_bfcp), 1),
    )

    # Retrieve the basis header information
    glp_get_bhead = cfunc(_glp+'get_bhead', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('k', c_int, 1),
    )

    # Retrieve row index in the basis header
    glp_get_row_bind = cfunc(_glp+'get_row_bind', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('i', c_int, 1),
    )

    # Retrieve column index in the basis header
    glp_get_col_bind = cfunc(_glp+'get_col_bind', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('j', c_int, 1),
    )

    # Perform forward transformation
    glp_ftran = cfunc(_glp+'ftran', None,
        ('lp', POINTER(glp_prob), 1),
        ('x', c_double_p, 1),
    )

    # Perform backward transformation
    glp_btran = cfunc(_glp+'btran', None,
        ('lp', POINTER(glp_prob), 1),
        ('x', c_double_p, 1),
    )



if _version >= (4, 9):
    # "warm up" LP basis
    lpx_warm_up = cfunc(_lpx+'warm_up', c_int,
        ('lp', POINTER(LPX), 1),
    )

if _version >= (4, 9):
    # compute row of the simplex table
    lpx_eval_tab_row = cfunc(_lpx+'eval_tab_row', c_int,
        ('lp', POINTER(LPX), 1),
        ('k', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # compute column of the simplex table
    lpx_eval_tab_col = cfunc(_lpx+'eval_tab_col', c_int,
        ('lp', POINTER(LPX), 1),
        ('k', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

if _version >= (4, 25):
    # compute row of the simplex tableau
    glp_eval_tab_row = cfunc(_glp+'eval_tab_row', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('k', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # compute column of the simplex tableau
    glp_eval_tab_col = cfunc(_glp+'eval_tab_col', c_int,
        ('lp', POINTER(glp_prob), 1),
        ('k', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

if _version >= (4, 9):
    # transform explicitly specified row
    lpx_transform_row = cfunc(_lpx+'transform_row', c_int,
        ('lp', POINTER(LPX), 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # transform explicitly specified column
    lpx_transform_col = cfunc(_lpx+'transform_col', c_int,
        ('lp', POINTER(LPX), 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
    )

    # perform primal ratio test
    lpx_prim_ratio_test = cfunc(_lpx+'prim_ratio_test', c_int,
        ('lp', POINTER(LPX), 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
        ('how', c_int, 1),
        ('tol', c_double, 1),
    )

    # perform dual ratio test
    lpx_dual_ratio_test = cfunc(_lpx+'dual_ratio_test', c_int,
        ('lp', POINTER(LPX), 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
        ('how', c_int, 1),
        ('tol', c_double, 1),
    )


#=============================================================================
# Branch-and-cut interface routines
#=============================================================================

if _version >= (4, 20):
    # reason codes:
    GLP_IROWGEN  = 0x01  # request for row generation
    GLP_IBINGO   = 0x02  # better integer solution found
    GLP_IHEUR    = 0x03  # request for heuristic solution
    GLP_ICUTGEN  = 0x04  # request for cut generation

if _version >= (4, 21):
    GLP_IBRANCH  = 0x05  # request for branching
    GLP_ISELECT  = 0x06  # request for subproblem selection
    GLP_IPREPRO  = 0x07  # request for preprocessing

if _version >= (4, 20):
    # Determine reason for calling the callback routine
    glp_ios_reason = cfunc(_glp+'ios_reason', c_int,
        ('tree', POINTER(glp_tree), 1),
    )

    # Access the problem object
    glp_ios_get_prob = cfunc(_glp+'ios_get_prob', POINTER(glp_prob),
        ('tree', POINTER(glp_tree), 1),
    )

    # Determine size of the branch-and-bound tree
    glp_ios_tree_size = cfunc(_glp+'ios_tree_size', None,
        ('tree', POINTER(glp_tree), 1),
        ('a_cnt', c_int_p, 1),
        ('n_cnt', c_int_p, 1),
        ('t_cnt', c_int_p, 1),
    )

    # Determine current active subproblem
    glp_ios_curr_node = cfunc(_glp+'ios_curr_node', c_int,
        ('tree', POINTER(glp_tree), 1),
    )

    # Determine next active subproblem
    glp_ios_next_node = cfunc(_glp+'ios_next_node', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

    # Determine previous active subproblem
    glp_ios_prev_node = cfunc(_glp+'ios_prev_node', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

    # Determine parent subproblem
    glp_ios_up_node = cfunc(_glp+'ios_up_node', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

    # Determine subproblem level
    glp_ios_node_level = cfunc(_glp+'ios_node_level', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

    # Determine subproblem local bound
    glp_ios_node_bound = cfunc(_glp+'ios_node_bound', c_double,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

    # Find active subproblem with best local bound
    glp_ios_best_node = cfunc(_glp+'ios_best_node', c_int,
        ('tree', POINTER(glp_tree), 1),
    )

    # Compute relative MIP gap
    glp_ios_mip_gap = cfunc(_glp+'ios_mip_gap', c_double,
        ('tree', POINTER(glp_tree), 1),
    )

if _version >= (4, 21):
    # Access subproblem application-specific data
    glp_ios_node_data = cfunc(_glp+'ios_node_data', c_void_p,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )

if _version >= (4, 32):
    def GLP_ATTR_FIELDS():
        return [ # additional row attributes
            ('level', c_int), # subproblem level at which the row was added
            ('origin', c_int), # the row origin flag
            ('klass', c_int), # the row class descriptor
            ('foo_bar', c_double*7), # (reserved)
        ]

    GLP_RF_REG      = 0  # regular constraint
    GLP_RF_LAZY     = 1  # "lazy" constraint
    GLP_RF_CUT      = 2  # cutting plane constraint

    GLP_RF_GMI      = 1  # Gomory's mixed integer cut
    GLP_RF_MIR      = 2  # mixed integer rounding cut
    GLP_RF_COV      = 3  # mixed cover cut
    GLP_RF_CLQ      = 4  # clique cut

    class glp_attr(Structure):
        _fields_ = GLP_ATTR_FIELDS()

    # retrieve additional row attributes
    glp_ios_row_attr = cfunc(_glp+'ios_row_attr', None,
        ('tree', POINTER(glp_tree), 1),
        ('i', c_int, 1),
        ('attr', POINTER(glp_attr), 1),
    )

    # determine current size of the cut pool
    glp_ios_pool_size = cfunc(_glp+'ios_pool_size', c_int,
        ('tree', POINTER(glp_tree), 1),
    )

    # add row (constraint) to the cut pool
    glp_ios_add_row = cfunc(_glp+'ios_add_row', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('name', c_char_p, 1),
        ('klass', c_int, 1),
        ('flags', c_int, 1),
        ('len', c_int, 1),
        ('ind', c_int_p, 1),
        ('val', c_double_p, 1),
        ('rhs', c_double, 1),
    )

    # remove row (constraint) from the cut pool
    glp_ios_del_row = cfunc(_glp+'ios_del_row', None,
        ('tree', POINTER(glp_tree), 1),
        ('i', c_int, 1),
    )

    # remove all rows (constraints) from the cut pool
    glp_ios_clear_pool = cfunc(_glp+'ios_clear_pool', None,
        ('tree', POINTER(glp_tree), 1),
    )

if _version >= (4, 21):
    # Check if can branch upon specified variable
    glp_ios_can_branch = cfunc(_glp+'ios_can_branch', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('j', c_int, 1),
    )

    # Choose variable to branch upon
    glp_ios_branch_upon = cfunc(_glp+'ios_branch_upon', None,
        ('tree', POINTER(glp_tree), 1),
        ('j', c_int, 1),
        ('sel', c_int, 1),
    )

if _version >= (4, 33):
    # branch selection indicator:
    GLP_NO_BRNCH       = 0  # select no branch
    GLP_DN_BRNCH       = 1  # select down-branch
    GLP_UP_BRNCH       = 2  # select up-branch

if _version >= (4, 21):
    # Select subproblem to continue the search
    glp_ios_select_node = cfunc(_glp+'ios_select_node', None,
        ('tree', POINTER(glp_tree), 1),
        ('p', c_int, 1),
    )


if _version >= (4, 20):
    # Provide solution found by heuristic
    glp_ios_heur_sol = cfunc(_glp+'ios_heur_sol', c_int,
        ('tree', POINTER(glp_tree), 1),
        ('x', c_double_p, 1),
    )

if _version >= (4, 20):
    # Terminate the solution process
    glp_ios_terminate = cfunc(_glp+'ios_terminate', None,
        ('tree', POINTER(glp_tree), 1),
    )


#=============================================================================
# Library environment routines
#=============================================================================

if _version >= (4, 9) and _version <= (4, 15):
    glp_version = lpx_version

if _version >= (4, 16):
    # Determine library version
    glp_version = cfunc(_glp+'version', c_char_p,
    )

if _version >= (4, 16) and _version <= (4, 27):
    class glp_ulong(Structure):
        _fields_ = [
            ('lo', c_uint),
            ('hi', c_uint),
        ]

if _version >= (4, 28):
    class glp_long(Structure):
        _fields_ = [
            ('lo', c_int),
            ('hi', c_int),
        ]

if _version >= (4, 21):
    # Enable/disable terminal output
    glp_term_out = cfunc(_glp+'term_out', None,
        ('flag', c_int, 1),
    )

if _version >= (4, 16):
    # Terminal hook function
    glp_term_hook_func = CFUNCTYPE(c_int, c_void_p, c_char_p)

    # Determine library version
    glp_term_hook = cfunc(_glp+'term_hook', None,
        ('func', glp_term_hook_func, 1),
        ('info', c_void_p, 1),
    )

if _version >= (4, 16) and _version <= (4, 27):
    # Get memory usage information
    glp_mem_usage = cfunc(_glp+'mem_usage', None,
        ('count', c_int_p, 1),
        ('cpeak', c_int_p, 1),
        ('total', POINTER(glp_ulong), 1),
        ('tpeak', POINTER(glp_ulong), 1),
    )

if _version >= (4, 28):
    # Get memory usage information
    glp_mem_usage = cfunc(_glp+'mem_usage', None,
        ('count', c_int_p, 1),
        ('cpeak', c_int_p, 1),
        ('total', POINTER(glp_long), 1),
        ('tpeak', POINTER(glp_long), 1),
    )

if _version >= (4, 19):
    # Set memory usage limit
    glp_mem_limit = cfunc(_glp+'mem_limit', None,
        ('limit', c_int, 1),
    )

if _version >= (4, 23):
    # Free GLPK library environment
    glp_free_env = cfunc(_glp+'free_env', None,
    )


#=============================================================================
# Wrap up all the functions and constants into __all__
#=============================================================================
__all__ = [x for x in locals().keys() \
    if  x.startswith('lpx_') or \
        x.startswith('LPX') or \
        x.startswith('glp') or \
        x.startswith('GLP')]

if __name__ == "__main__":
    print "Welcome. You are using ctypes-glpk, a Python wrapper for GLPK written by Minh-Tri Pham."
