"""
CBMPy: CBWrite module
=====================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2024 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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

Author: Brett G. Olivier PhD
Contact developers: https://github.com/SystemsBioinformatics/cbmpy/issues
Last edit: $Author: bgoli $ ($Id: CBWrite.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals

if 'cDir' in vars():
    cDir = vars()['cDir']
else:
    cDir = None
import os, time, numpy, zipfile

# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio

from . import CBTools, CBXML


from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

from packaging import version as pkgver

_HAVE_SYMPY_ = False
try:
    import sympy

    if pkgver.parse(sympy.__version__) >= pkgver.Version('0.7.5'):
        _HAVE_SYMPY_ = True

    else:
        del sympy
        print(
            '\nWARNING: SymPy version 0.7.5 or newer is required for symbolic matrix support.'
        )
except ImportError:
    _HAVE_SYMPY_ = False


_HAVE_XLWT_ = False
try:
    import xlwt

    _HAVE_XLWT_ = True
except ImportError:
    print('\nINFO: No xlwt module available, Excel spreadsheet creation disabled')


def saveModel(model, filename, compress=False):
    """
    Saves the model to an SBML file using the latest SBML3 FBC version.

    Parameters
    ----------
    model : CBMPy model
        The CBMPy model to be saved.
    filename : str
        The destination filename for the SBML file.
    compress : bool, optional
        If True, compresses the model using ZIP encoding. Default is False.

    Example
    --------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('path_to_model.xml')
    >>> CBWrite.saveModel(model, 'saved_model.xml', compress=True)
    """

    writeSBML3FBCV2(
        model,
        filename,
        directory=None,
        gpr_from_annot=False,
        add_groups=True,
        add_cbmpy_annot=True,
        add_cobra_annot=False,
        validate=False,
        compress_bounds=True,
        zip_model=compress,
        return_model_string=False,
    )


def writeSBML3FBC(
    fba,
    fname,
    directory=None,
    gpr_from_annot=False,
    add_groups=True,
    add_cbmpy_annot=True,
    add_cobra_annot=False,
    xoptions={'fbc_version': 1, 'validate': False, 'compress_bounds': True},
):
    """
    Takes an FBA model object and writes it to file as SBML L3 FBC.

    Parameters
    ----------
    fba : CBMPy model object
        The FBA model to be written to file.
    fname : str
        The filename for the output SBML file.
    directory : str, optional
        The directory where the output file will be written. If not specified, uses current directory.
    gpr_from_annot : bool, optional
        If True, will attempt to add gene-protein associations from model annotations.
    add_groups : bool, optional
        If True, adds SBML Level 3 groups support if available in libSBML.
    add_cbmpy_annot : bool, optional
        If True, adds CBMPy-specific KeyValueData annotations.
    add_cobra_annot : bool, optional
        If True, adds COBRA-style notes from the annotations.
    xoptions : dict, optional
        A dictionary of extended options including keys for 'fbc_version', 'validate',
        'compress_bounds', 'zip_model', and 'return_model_string'.

    Returns
    -------
    None or str
        None by default. If 'return_model_string' in xoptions is True, returns the SBML string.

    Example
    -------
    >>> from cbmpy import CBModel, CBWrite
    >>> model = CBModel.load_model('example_model.xml')
    >>> CBWrite.writeSBML3FBC(model, 'output_model.xml', directory='/path/to/output', xoptions={'fbc_version': 2, 'validate': True})
    """
    sbml_level_version = (3, 1)
    autofix = True
    return_fbc = False
    return CBXML.sbml_writeSBML3FBC(
        fba,
        fname,
        directory,
        sbml_level_version,
        autofix,
        return_fbc,
        gpr_from_annot,
        add_groups,
        add_cbmpy_annot,
        add_cobra_annot,
        xoptions,
    )


def writeSBML3FBCV2(
    fba,
    fname,
    directory=None,
    gpr_from_annot=False,
    add_groups=True,
    add_cbmpy_annot=True,
    add_cobra_annot=False,
    validate=False,
    compress_bounds=False,
    zip_model=False,
    return_model_string=False,
):
    """
    Writes an FBA (Flux Balance Analysis) model object to a file in SBML (Systems Biology Markup Language) Level 3 Version 1 Format with FBC (Flux Balance Constraints) Version 2.

    Parameters
    ----------
    fba : object
        The FBA model object to write to file.
    fname : str
        The file name to write the SBML to.
    directory : str, optional
        The directory where to save the file. Defaults to the current directory if None.
    gpr_from_annot : bool, optional
        If True, gene-protein associations will be added from annotations.
    add_groups : bool, optional
        If True, SBML level 3 groups will be added if supported by libSBML.
    add_cbmpy_annot : bool, optional
        If True, CBMPy KeyValueData annotations will be added.
    add_cobra_annot : bool, optional
        If True, COBRA-style annotations will be added.
    validate : bool, optional
        If True, the generated SBML file will be validated.
    compress_bounds : bool, optional
        If True, attempts will be made to compress flux bound parameters.
    zip_model : bool, optional
        If True, the model will be compressed using ZIP encoding.
    return_model_string : bool, optional
        If True, returns the SBML file as a string instead of writing to disk.

    Returns
    -------
    None or str
        None by default; if 'return_model_string' is True, returns the generated SBML as a string.

    Example
    -------
    >>> fba_model = CBModel.loadModel('example_model.xml')
    >>> writeSBML3FBCV2(fba_model, 'example_model_fbcv2.xml', directory='./models/', add_cbmpy_annot=True, add_cobra_annot=False, validate=True, zip_model=False)
    """

    xoptions = {
        'fbc_version': 2,
        'validate': validate,
        'compress_bounds': compress_bounds,
        'return_model_string': return_model_string,
        'zip_model': zip_model,
    }
    sbml_level_version = (3, 1)
    autofix = (True,)
    return_fbc = False
    return CBXML.sbml_writeSBML3FBC(
        fba,
        fname,
        directory,
        sbml_level_version,
        autofix,
        return_fbc,
        gpr_from_annot,
        add_groups,
        add_cbmpy_annot,
        add_cobra_annot,
        xoptions,
    )


def writeSBML3FBCV3(
    fba,
    fname,
    directory=None,
    gpr_from_annot=False,
    add_groups=True,
    add_cbmpy_annot=True,
    add_cobra_annot=False,
    validate=False,
    compress_bounds=False,
    zip_model=False,
    return_model_string=False,
):
    """
    Writes an FBA model object to a file in SBML Level 3 Format with FBC version 3.

    Parameters
    ----------
    fba : CBMPy model object
        The FBA model to be written to file.
    fname : str
        The filename for the output SBML file.
    directory : str, optional
        The directory where the output file will be written. If not specified, uses current directory.
    gpr_from_annot : bool, optional
        If True, will attempt to add gene-protein associations from model annotations.
    add_groups : bool, optional
        If True, adds SBML Level 3 groups support if available in libSBML.
    add_cbmpy_annot : bool, optional
        If True, adds CBMPy-specific KeyValueData annotations.
    add_cobra_annot : bool, optional
        If True, adds COBRA-style notes from the annotations.
    validate : bool, optional
        If True, the generated SBML file will be validated.
    compress_bounds : bool, optional
        If True, attempts will be made to compress flux bound parameters.
    zip_model : bool, optional
        If True, the model will be compressed using ZIP encoding.
    return_model_string : bool, optional
        If True, returns the SBML file as a string instead of writing to disk.

    Returns
    -------
    None or str
        None by default; if 'return_model_string' is True, returns the generated SBML as a string.

    Example
    -------
    >>> from cbmpy import CBModel, CBWrite
    >>> model = CBModel.load_model('example_model.xml')
    >>> CBWrite.writeSBML3FBCV3(model, 'output_model_sbmlv3.xml', directory='/path/to/output', add_cbmpy_annot=True, validate=True)
    """

    xoptions = {
        'fbc_version': 3,
        'validate': validate,
        'compress_bounds': compress_bounds,
        'return_model_string': return_model_string,
        'zip_model': zip_model,
    }
    sbml_level_version = (3, 1)
    autofix = (True,)
    return_fbc = False
    if xoptions['fbc_version'] == 3:
        add_cobra_annot = False
    return CBXML.sbml_writeSBML3FBC(
        fba,
        fname,
        directory,
        sbml_level_version,
        autofix,
        return_fbc,
        gpr_from_annot,
        add_groups,
        add_cbmpy_annot,
        add_cobra_annot,
        xoptions,
    )


def writeCOBRASBML(fba, fname, directory=None):
    """
    Write an FBA model object to a file in a pre-2010 COBRA-compatible dialect.

    Parameters
    ----------
    fba : CBMPy model
        The FBA model object to write to file.
    fname : str
        The filename of the file to write to.
    directory : str, optional
        The directory where to save the file. If not specified, uses current directory.

    Examples
    --------
    >>> fba_model = CBMPy.readSBML3FBC('my_model.xml')
    >>> CBWrite.writeCOBRASBML(fba_model, 'cobra_model.xml', directory='./models')

    """
    return CBXML.sbml_writeCOBRASBML(fba, fname, directory)


def writeSBML2FBA(fba, fname, directory=None, sbml_level_version=None):
    """
    Writes an FBA model object to an SBML L2 file with FBA annotations.

    Parameters
    ----------
    fba : CBMPy model object
        The FBA model to be written to the SBML file.
    fname : str
        The name of the file where the SBML model will be written.
    directory : str, optional
        The directory where the file will be saved. If not specified, the file
        will be saved in the current working directory.
    sbml_level_version : tuple of int, optional
        Specifies the SBML level and version. For BiGG/FAME style annotations,
        use (2, 1). If not specified, the default SBML level and version will be used.

    Example
    -------
    >>> fba_model = CBMPy.readSBML3FBC('path_to\model.xml')
    >>> writeSBML2FBA(fba_model, 'new_model.xml', directory='path_to\save_directory', sbml_level_version=(2, 1))

    """
    CBXML.sbml_writeSBML2FBA(fba, fname, directory, sbml_level_version)


def writeSensitivitiesToCSV(sensitivities, fname):
    """
    Writes out a sensitivity report to CSV files, detailing objective sensitivities,
    and bound sensitivities. The output comprises two files:
    one for flux sensitivities and another for constraint sensitivities.

    Parameters
    ----------
    sensitivities : tuple
        A tuple containing dictionaries of objective coefficient sensitivities
        (`obj_sens`), constraint RHS sensitivities (`rhs_sens`), and
        bound sensitivities (`bound_sens`).
    fname : str
        Base filename for the output CSV files. The function generates
        two files appending `_flux_sensitivity.csv` and `_constraint_sensitivity.csv`
        to this base name.

    Example
    -------
    >>> sensitivities = (obj_sens, rhs_sens, bound_sens)
    >>> writeSensitivitiesToCSV(sensitivities, 'model_sensitivities')
    This will create 'model_sensitivities_flux_sensitivity.csv' and
    'model_sensitivities_constraint_sensitivity.csv'.

    """
    obj_sens = sensitivities[0]
    rhs_sens = sensitivities[1]
    bound_sens = sensitivities[2]
    F = open(fname + '_flux_sensitivity.csv', 'w')
    head = "Flux,Reduced cost,OCS low,OC value,OCS high,LB low,LB high,UB low,UB high"
    F.write(head + '\n')
    for j in obj_sens:
        rc = obj_sens[j][0]
        lcs = obj_sens[j][1]
        ocv = obj_sens[j][2]
        ucs = obj_sens[j][3]
        if j in bound_sens:
            lbs = bound_sens[j][0]
            lb = bound_sens[j][1]
            ub = bound_sens[j][2]
            ubs = bound_sens[j][3]
        else:
            lbs = 0
            lb = 0
            ub = 0
            ubs = 0
        F.write(
            '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (j, rc, lcs, ocv, ucs, lbs, lb, ub, ubs)
        )

    for j in bound_sens:
        if j not in obj_sens:
            rc = 0
            lcs = 0
            ocv = 0
            ucs = 0
            lbs = bound_sens[j][0]
            lb = bound_sens[j][1]
            ub = bound_sens[j][2]
            ubs = bound_sens[j][3]
            F.write(
                '%s,%s,%s,%s,%s,%s,%s,%s,%s\n'
                % (j, rc, lcs, ocv, ucs, lbs, lb, ub, ubs)
            )
    F.flush()
    F.close()

    F = open(fname + '_constraint_sensitivity.csv', 'w')
    F.write('Constraint,RHS low,RHS,RHS high\n')
    for c in rhs_sens:
        F.write('%s,%s,%s,%s\n' % (c, rhs_sens[c][0], rhs_sens[c][1], rhs_sens[c][2]))
    F.flush()
    F.close()


def WriteModelRaw(fba, work_dir=cDir):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelRaw()\"

    This function acts as a wrapper for the `writeModelRaw()` function, alerting the user about the deprecation
    of this naming convention and suggesting the updated function name.

    Parameters
    ----------
    fba : CBMPy model
        An instantiated FBAmodel instance.
    work_dir : str, optional
        Directory designated for output. If None, uses the current directory.

    Example
    -------
    >>> from cbmpy import CBModel, CBWrite
    >>> my_model = CBModel.CBRead.readSBML3FBC('path_to_my_model.xml')
    >>> CBWrite.WriteModelRaw(my_model, work_dir='./model_output/')
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeModelRaw()\"\n'
    )
    time.sleep(1)
    writeModelRaw(fba, work_dir)


def writeModelRaw(fba, work_dir=cDir):
    """
    Write the model data to a raw text format for inspection or debugging.

    This function provides a simple textual dump of the basic components of a CBMPy model including species, reactions, flux bounds, objectives, and stoichiometry.

    Parameters
    ----------
    fba : CBMPy.FBA
        The CBMPy FBA model instance to be dumped.
    work_dir : str, optional
        The working directory path where the output file will be saved. Defaults to the current working directory.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('example_model.xml')
    >>> CBWrite.writeModelRaw(model, work_dir='output_directory')

    """
    if work_dir == None:
        work_dir = os.getcwd()
    FF = open(os.path.join(work_dir, 'WriteModelRawOutput.txt'), 'w')
    FF.write('Species information:\n\n')
    for s in fba.species:
        FF.write(
            '%s: value=%f, is_boundary=%s, name=%s\n'
            % (s.getId(), s.value, int(s.is_boundary), s.name)
        )
    FF.write('\nReaction information:\n\n')
    for r in fba.reactions:
        FF.write(
            '%s: reversible=%s, name=%s\n' % (r.getId(), int(r.reversible), r.name)
        )
        FF.write('\t%s\n' % r.getSpeciesIds())
    FF.write('\nConstraint information:\n\n')
    for c in fba.flux_bounds:
        FF.write('%s %s %f\n' % (c.reaction, c.operation, c.value))
    FF.write('\nObjective information:\n\n')
    FF.write(
        'Active Objective: %s (%s)\n'
        % (
            fba.objectives[fba.activeObjIdx].getId(),
            fba.objectives[fba.activeObjIdx].operation,
        )
    )
    for o in fba.objectives:
        FF.write(
            '%s: %s, %s\n' % (o.getId(), o.operation, o.getFluxObjectiveReactions())
        )
    if hasattr(fba.N, 'shape'):
        FF.write('\nStoichiometric information:\n\n')
        FF.write('N-matrix dimensions = (%s,%s)\n' % fba.N.shape)
    FF.close()
    print(
        'WriteModelRaw has written a file to {}'.format(
            os.path.join(work_dir, 'WriteModelRawOutput.txt')
        )
    )


def BuildLPFluxBounds(fba, use_rational=False):
    """
    Generates a string buffer containing flux bounds for each reaction in the FBA model formatted for linear programming (LP) constraints.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model whose flux bounds are to be formatted for LP.
    use_rational : bool, optional
        If True, exports bounds as rational numbers. Requires SymPy. Default is False.

    Returns
    -------
    cStringIO.StringIO
        A buffer containing the LP-formatted flux bounds for each reaction.

    Notes
    -----
    This function creates flux bounds in the form "reactionID <= value" for upper bounds or "reactionID >= value" for lower bounds. It is designed to integrate with LP solvers by providing a mechanism to incorporate flux bounds directly.

    Example
    -------
    >>> fba_model = CBMPy.readSBML3FBC('example_model.xml')
    >>> flux_bounds_buffer = BuildLPFluxBounds(fba_model)
    >>> print(flux_bounds_buffer.getvalue())
    """
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    FFS = csio.StringIO()
    c2s = {}
    for c in fba.flux_bounds:
        if __DEBUG__:
            print('%s: %s' % (c.getId(), c.operation))
        minv = None
        maxv = None
        oper = None
        R = c.reaction
        if c.operation in ['less', 'lessEqual']:
            maxv = c.value
            oper = c.operation
            if oper == 'less':
                oper = '<='
                if __DEBUG__:
                    print('LP Bounds operator must be <= not <')
            else:
                oper = '<='
        if c.operation in ['greater', 'greaterEqual']:
            minv = c.value
            oper = c.operation
            if oper == 'greater':
                oper = '<='
                if __DEBUG__:
                    print('LP Bounds operator must be <= not <')
            else:
                oper = '<='
        if c.operation == 'equals':
            maxv = c.value
            oper = '='
        if maxv != None:
            if use_rational:
                c2s.update({R + 'r': '%s %s' % (oper, sympy.Rational(maxv))})
            else:
                if maxv == float('inf'):
                    c2s.update({R + 'r': '%s +%s' % (oper, maxv)})
                else:
                    c2s.update({R + 'r': '%s %s' % (oper, maxv)})
        elif minv != None:
            if use_rational:
                c2s.update({R + 'l': '%s %s' % (sympy.Rational(minv), oper)})
            else:
                c2s.update({R + 'l': '%s %s' % (minv, oper)})

    c2sk = list(c2s)
    c2sk.sort()
    for r in c2sk:
        rht = None
        lht = None
        R = r[:-1]
        if R + 'r' in c2s:
            rht = c2s[R + 'r']
            c2s.pop(R + 'r')
        if R + 'l' in c2s:
            lht = c2s[R + 'l']
            c2s.pop(R + 'l')
        if lht != None and rht != None:
            if __DEBUG__:
                print('{} {} {}'.format(lht, R, rht))
            ##  FFS.write('%s: %s %s %s\n' % (R, lht, R, rht))
            FFS.write('%s %s %s\n' % (lht, R, rht))  # THIS MUST BE SO for GLPK
        elif lht == None and rht == None:
            if __DEBUG__:
                print('Skipping: ({}, {}, {})'.format(lht, R, rht))
        elif lht != None:
            if __DEBUG__:
                print('{} {}'.format(lht, R))
            # FFS.write('%s: %s %s\n' % (R, lht, R))
            FFS.write('%s %s\n' % (lht, R))
        elif rht != None:
            if __DEBUG__:
                print('{} {}'.format(R, rht))
            # FFS.write('%s: %s %s\n' % (R, R, rht))
            FFS.write('%s %s\n' % (R, rht))
        else:
            print('CONFUSION: ({}, {}, {})'.format(lht, R, rht))
    return FFS


def BuildLPConstraints(fba, use_rational=False):
    """
    Constructs a StringIO object that contains the linear programming constraints
    derived from the FBA model's stoichiometry matrix.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model which has a stoichiometry matrix defined.
    use_rational : bool, optional
        If True, the constraints will be output as rational numbers. This requires
        SymPy to be installed. The default is False, which outputs floating-point numbers.

    Returns
    -------
    csio.StringIO
        A StringIO object that contains the model's constraints in a format suitable
        for linear programming solvers. Each constraint is written on a separate line.

    Examples
    --------
    >>> fba_model = CBMPy.readSBML3FBC('example_model.xml')
    >>> constraints_buffer = BuildLPConstraints(fba_model, use_rational=True)
    >>> print(constraints_buffer.getvalue())

    Notes
    -----
    If the model's stoichiometry matrix is inconsistent or undefined, the function attempts
    to rebuild it using `CBTools.addStoichToFBAModel(fba)`. If `use_rational` is True but SymPy
    is not installed, a warning will be issued and floating-point numbers will be used instead.
    """

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install SymPy for rational output.')
    rebuild_stoich = False
    if not hasattr(fba, 'N') or fba.N == None:
        rebuild_stoich = True
    else:
        if len([s for s in fba.species if not s.is_boundary]) != fba.N.array.shape[0]:
            rebuild_stoich = True
        elif len(fba.reactions) != fba.N.array.shape[1]:
            rebuild_stoich = True
    if rebuild_stoich:
        print(
            '\nWarning: FBA object has inconsistent stoichiometric matrix, rebuilding it now.'
        )
        CBTools.addStoichToFBAModel(fba)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r, c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName: RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in range(len(fba.N.row)):
        if len(constr[fba.N.row[r]]) > 0:
            FFS.write(' %s: ' % fba.N.row[r])
            for col in constr[fba.N.row[r]]:
                if use_rational:
                    if col[0] > 0.0:
                        FFS.write('+%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                    else:
                        FFS.write('%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                else:
                    if col[0] > 0.0:
                        if col[0] == 1.0:
                            FFS.write('+ %s ' % (col[1]))
                        else:
                            FFS.write('+%.20f %s ' % (col[0], col[1]))
                    else:
                        if col[0] == -1.0:
                            FFS.write('- %s ' % (col[1]))
                        else:
                            FFS.write('%.20f %s ' % (col[0], col[1]))
            operator = fba.N.operators[r].upper()
            if operator == 'E' or operator == '=':
                operator = '='
            elif operator == 'G' or operator == '>' or operator == '>=':
                operator = '>='
            elif operator == 'L' or operator == '<' or operator == '<=':
                operator = '<='
            else:
                raise RuntimeError('\nINVALID operator: %s' % operator)
            FFS.write('%s ' % operator)
            FFS.write('%s\n' % fba.N.RHS[r])
    return FFS


def BuildLPUserConstraints(fba, use_rational=False):
    """
    Construct a StringIO object that represents the user-defined constraints for an FBA model in linear programming format.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model for which user-defined constraints are constructed.
    use_rational : bool, optional
        If True, coefficients are represented as rational numbers. Requires SymPy. Defaults to False, using floating-point representation.

    Returns
    -------
    csio.StringIO
        A StringIO object that lists each user-defined constraint in linear programming format.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> constraint_buffer = CBWrite.BuildLPUserConstraints(fba_model, use_rational=True)
    >>> print(constraint_buffer.getvalue())

    Notes
    -----
    This function assumes your model has user-defined constraints (e.g., additional linear equations)
    set up in its 'CM' attribute, usually defined through CBMPy's constraint management API.

    """

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational output')
    rebuild_stoich = False
    assert fba.user_constraints != None, "\nNo user constraints to build"
    if not hasattr(fba, 'CM') or fba.CM == None:
        rebuild_stoich = True

    if rebuild_stoich:
        print(
            '\nWarning: FBA object has inconsistent user constraint matrix, rebuilding it now.'
        )
        CBTools.addStoichToFBAModel(fba)

    constr = {}
    for r in range(fba.CM.array.shape[0]):
        rowName = fba.CM.row[r]
        RCon = []
        for c in range(fba.CM.array.shape[1]):
            colName = fba.CM.col[c]
            colCoef = fba.CM.array[r, c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName: RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in range(len(fba.CM.row)):
        if len(constr[fba.CM.row[r]]) > 0:
            FFS.write(' %s: ' % fba.CM.row[r])
            for col in constr[fba.CM.row[r]]:
                if use_rational:
                    if col[0] > 0.0:
                        FFS.write('+%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                    else:
                        FFS.write('%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                else:
                    if col[0] > 0.0:
                        if col[0] == 1.0:
                            FFS.write('+ %s ' % (col[1]))
                        else:
                            FFS.write('+%.20f %s ' % (col[0], col[1]))
                    else:
                        if col[0] == -1.0:
                            FFS.write('- %s ' % (col[1]))
                        else:
                            FFS.write('%.20f %s ' % (col[0], col[1]))
            operator = fba.CM.operators[r]
            if operator == 'E' or operator == '=':
                operator = '='
            elif operator == 'G' or operator == '>' or operator == '>=':
                operator = '>='
            elif operator == 'L' or operator == '<' or operator == '<=':
                operator = '<='
            else:
                raise RuntimeError('\nINVALID operator: %s' % operator)
            FFS.write('%s ' % operator)
            FFS.write('%s\n' % fba.CM.RHS[r])
    return FFS


def BuildLPConstraintsRelaxed(fba):
    """
    Builds and returns a StringIO object containing the constraints for Flux Balance Analysis (FBA) in
    LP format where the constraints are relaxed to consider steady-state conditions as inequalities (dS/dt >= 0).

    Parameters
    ----------
    fba : CBMPy FBA model object
        The constraint-based model for which the relaxed LP constraints are to be generated.

    Returns
    -------
    csio.StringIO
        A StringIO object containing the LP constraints where steady-state conditions are described using
        inequalities (>= 0) to indicate fluxes can only be non-negative.

    Raises
    ------
    DeprecationWarning
        If this function is called, indicating it is deprecated.

    Example
    --------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> constraints_buffer = CBWrite.BuildLPConstraintsRelaxed(fba_model)
    >>> print(constraints_buffer.getvalue())

    Notes
    -----
    This function is deprecated and may be removed in future releases. It is recommended to use
    alternative methods for generating LP constraints.

    """
    raise DeprecationWarning("\nThis method is deprecated")
    if not hasattr(fba, 'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r, c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName: RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write('%s: ' % r)
            for col in constr[r]:
                if col[0] > 0.0:
                    FFS.write('+%.20f %s ' % (col[0], col[1]))
                else:
                    FFS.write('%.20f %s ' % (col[0], col[1]))
            FFS.write('>= 0\n')

    return FFS


def BuildLPConstraintsStrict(fba, use_rational=False):
    """
    Build and return a csio.StringIO object containing the constraints in LP format for strict
    stoichiometric constraints, where stoichiometric balances are assumed to be exactly zero (dS/dt = 0).

    Parameters
    ----------
    fba : FBA model object
        The FBA model for which the LP constraints are generated.
    use_rational : bool, optional
        If true, generates the constraints using rational numbers. Requires sympy. Default is False.

    Returns
    -------
    csio.StringIO
        A StringIO object containing the generated constraints in LP format.

    Raises
    ------
    RuntimeError
        If sympy is required (use_rational=True) but not installed.

    Examples
    --------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> constraints_buffer = CBWrite.BuildLPConstraintsStrict(fba_model, use_rational=True)
    >>> print(constraints_buffer.getvalue())

    Notes
    -----
    The output is compatible with various LP solvers that accept strict stoichiometric constraints.
    Ensure the FBA model has a stoichiometry matrix (N) defined, or has been pre-processed using CBTools.

    """
    ##  print "Consider using the new BuildLPConstraints() method"

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    if not hasattr(fba, 'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r, c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName: RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write(' %s: ' % r)
            for col in constr[r]:
                if use_rational:
                    if col[0] > 0.0:
                        FFS.write('+%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                    else:
                        FFS.write('%s %s ' % (sympy.Rational('%s' % col[0]), col[1]))
                else:
                    if col[0] > 0.0:
                        if col[0] == 1.0:
                            FFS.write('+ %s ' % (col[1]))
                        else:
                            FFS.write('+%.20f %s ' % (col[0], col[1]))
                    else:
                        if col[0] == -1.0:
                            FFS.write('- %s ' % (col[1]))
                        else:
                            FFS.write('%.20f %s ' % (col[0], col[1]))
            if len(constr[r]) != 0:
                FFS.write('= 0\n')
            else:
                FFS.write('\n')

    return FFS


def BuildLPConstraintsMath(fba, use_rational=False):
    """
    Construct a StringIO object that encodes the model constraints in LP format considering both directions
    where dS/dt => 0 (n1) and dS/dt <= 0 (n2) for mathematical programming formulations.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model for which the constraints are to be generated.
    use_rational : bool, optional
        If true, generates the constraints using rational numbers if the sympy library is available.
        Defaults to False, which uses floating-point numbers.

    Returns
    -------
    csio.StringIO
        A StringIO object that lists the generated constraints in LP format suitable for linear programming solvers.

    Raises
    ------
    DeprecationWarning
        Alerts that this particular function is deprecated.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> constraints_buffer = CBWrite.BuildLPConstraintsMath(fba_model, use_rational=True)
    >>> print(constraints_buffer.getvalue())

    Notes
    -----
    This function is deprecated and may not accurately represent the latest format or standards for LP model constraints.
    """
    raise DeprecationWarning("\nThis method is deprecated")
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    if not hasattr(fba, 'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r, c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName: RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write('%sn1: ' % r)
            ##  FFS.write(' ')
            for col in constr[r]:
                if col[0] > 0.0:
                    FFS.write('+%s %s ' % (col[0], col[1]))
                else:
                    FFS.write('%s %s ' % (col[0], col[1]))
            # This is a fudge
            FFS.write('>= 0\n')

    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write('%sn2: ' % r)
            for col in constr[r]:
                if -col[0] > 0.0:
                    FFS.write('+%s %s ' % (-col[0], col[1]))
                else:
                    FFS.write('%s %s ' % (-col[0], col[1]))
            # This is a fudge
            FFS.write('>= 0\n')
    return FFS


def WriteModelLPOld(
    fba,
    work_dir=None,
    multisymb=' ',
    lpt=True,
    constraint_mode='strict',
    use_rational=False,
    format='%s',
):
    """
    A deprecated method for writing an FBA model to a file in LP format.

    This function acts as a wrapper for ``writeModelLPOld()``, signaling the function's deprecation.

    Parameters
    ----------
    fba : FBA model instance
        The FBA model to be written.
    work_dir : string, optional
        The directory where the output file will be saved. By default, uses the current directory.
    multisymb : string, optional
        The multiplication symbol used in the output. By default, it's a whitespace.
    lpt : bool, optional
        If True, the file will be saved in LPT format. By default, it's True.
    constraint_mode : string, optional
        The mode of constraint to be applied. By default, it's 'strict'.
    use_rational : bool, optional
        If True, rational numbers will be used. By default, it's False.
    format : string, optional
        The format string for numbers. By default, it's '%s'.

    Notes
    -----
    This function is deprecated and will be removed in future releases. Use ``writeModelLPOld()`` instead.

    Example
    -------
    >>> fba_model = cbmpy.CBRead.readSBML3FBC('model_file.xml')
    >>> WriteModelLPOld(fba_model, work_dir='output/', lpt=True, use_rational=True)
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeModelLPOld()\"\n'
    )
    time.sleep(1)
    writeModelLPOld(
        fba, work_dir, multisymb, lpt, constraint_mode, use_rational, format
    )


def writeModelLPOld(
    fba,
    work_dir=None,
    multisymb=' ',
    lpt=True,
    constraint_mode='strict',
    use_rational=False,
    format='%s',
):
    """
    Writes an FBA model as an LP or LPT format.

    Parameters
    ----------
    fba : FBA model instance
        The FBA model to be written to a file.
    work_dir : str, optional
        The directory where the output file will be saved. If None, the current directory is used.
    multisymb : str, optional
        The multiplication symbol used in formatting the LP file. Defaults to ' ' (space).
    lpt : bool, optional
        If True, generates an LPT file. If False, generates an LP file. Defaults to True.
    constraint_mode : {'strict', 'relaxed', 'math'}, optional
        The mode for generating constraints. Defaults to 'strict'.
    use_rational : bool, optional
        If True, numbers will be represented as rationals. This requires sympy to be installed.
        Defaults to False.
    format : str, optional
        The format specifier for floating-point numbers. Defaults to '%s'.

    Notes
    -----
    This function is deprecated. Consider using updated methods for LP file generation.

    """

    print("\nTHIS FUNCTION IS DEPRECATED\n")

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    FNAME = None
    if work_dir == None:
        FnameTmp = fba.getId()
    else:
        FnameTmp = os.path.join(work_dir, fba.getId())
    if use_rational:
        FnameTmp = FnameTmp + '.rat'
    if not lpt:
        FNAME = FnameTmp + '.lp'
        FF = open(FNAME, 'w')
        FF.write('Problem\n %s\n\n' % FnameTmp)
    else:
        FNAME = FnameTmp + '.lp'
        FF = open(FNAME, 'w')
        FF.write('\\\\ %s \n\n' % FnameTmp)
    objO = fba.objectives[fba.activeObjIdx].operation.lower()
    objO = objO[0].upper() + objO[1:]
    FF.write('%s\n' % objO)
    objStr = '%s_objf: ' % fba.objectives[fba.activeObjIdx].getId()
    for fObj in fba.objectives[fba.activeObjIdx].flux_objectives:
        sign = None
        nc = 0.0
        try:
            nc = float(fObj.coefficient)
        except ValueError:
            print(
                'Suspected rational number ({}) detected in fluxObjective {}'.format(
                    fObj.coefficient, fObj.getId()
                )
            )
        if nc >= 0.0:
            sign = '+'
        else:
            sign = '-'
        if use_rational:
            objStr += ' %s%s%s' % (
                sympy.Rational(fObj.coefficient),
                multisymb,
                fObj.reaction,
            )
        else:
            objStr += ' %s %s%s%s' % (sign, abs(nc), multisymb, fObj.reaction)
    FF.write('%s\n' % objStr)
    if constraint_mode == 'math':
        CONST = BuildLPConstraintsMath(fba)
    elif constraint_mode == 'relaxed':
        CONST = BuildLPConstraintsRelaxed(fba)
    else:
        CONST = BuildLPConstraintsStrict(fba, use_rational)

    CONST.seek(0)
    BOUNDS = BuildLPFluxBounds(fba, use_rational)
    BOUNDS.seek(0)
    if __DEBUG__:
        print(CONST.read())
        CONST.seek(0)
    if __DEBUG__:
        print(BOUNDS.read())
        BOUNDS.seek(0)

    FF.write('\nSubject To\n')
    FF.write(CONST.read())
    FF.write('\nBounds\n')
    FF.write(BOUNDS.read())
    FF.write('\nEND\n')
    FF.close()
    print('writeModelLP has written a file to {}'.format(FNAME))
    return FNAME


def WriteModelLP(
    fba,
    work_dir=None,
    fname=None,
    multisymb=' ',
    format='%s',
    use_rational=False,
    constraint_mode=None,
    quiet=False,
):
    """
    Writes the given FBA model to a file in LP format. This function is a legacy bridge
    for compatibility and directs to the updated method `writeModelLP` with updated arguments.

    Parameters
    ----------
    fba : CBMPy FBA model instance
        The FBA model to write.
    work_dir : str, optional
        The directory where to write the LP file. Defaults to the current directory.
    fname : str, optional
        The name of the file to write. If not specified, uses the model ID.
    multisymb : str, optional
        The symbol used to denote multiplication in the LP file. Defaults to ' '.
    format : str, optional
        The format string for numbers. Defaults to '%s'.
    use_rational : bool, optional
        If True, numbers are written as rationals. Requires sympy. Defaults to False.
    constraint_mode : str, optional
        This parameter is deprecated and has no effect.
    quiet : bool, optional
        If True, suppresses print statements. Defaults to False.

    Example
    --------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('example_model.xml')
    >>> CBWrite.WriteModelLP(model, fname='model.lp", work_dir='/path/to/dir')

    Notes
    -----
    This function wraps around the new `writeModelLP` to provide backward compatibility
    with older scripts that may rely on its name or functionality.
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeModelLP()\"\n'
    )
    time.sleep(1)
    writeModelLP(
        fba, work_dir, fname, multisymb, format, use_rational, constraint_mode, quiet
    )


def writeModelLP(
    fba,
    work_dir=None,
    fname=None,
    multisymb=' ',
    format='%s',
    use_rational=False,
    constraint_mode=None,
    quiet=False,
):
    """
    Writes an FBA object as an LP in CPLEX LP format.

    Parameters
    ----------
    fba : FBA model instance
        The FBA model to be written.
    work_dir : str, optional
        The directory where the output file will be saved. If None, the current directory is used.
    fname : str, optional
        The name of the file to write. If not specified, uses the model ID.
    multisymb : str, optional
        The multiplication symbol used in the LP file. Defaults to ' ' (space).
    format : str, optional
        The format string for numbers. Defaults to '%s'.
    use_rational : bool, optional
        If True, numbers are written as rationals. Requires sympy. Defaults to False.
    constraint_mode : str, optional
        This parameter is deprecated and has no effect.
    quiet : bool, optional
        If True, suppresses print statements. Defaults to False.

    Returns
    -------
    str
        The path to the written LP file.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('example_model.xml')
    >>> lp_file_path = CBWrite.writeModelLP(model, work_dir='/path/to/dir', quiet=True)
    >>> print(lp_file_path)
    """
    if constraint_mode != None:
        print("\nConstraint_mode has been deprecated")
        time.sleep(5)

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print(
            '\nWarning switching to floating point arithmetic: install Sympy for rational IO'
        )
        time.sleep(2)

    FNAME = None
    if fname == None:
        fname = fba.getId()
    if work_dir != None:
        FnameTmp = os.path.join(work_dir, fname)
    else:
        FnameTmp = fname

    FNAME = FnameTmp + '.lp'
    FF = open(FNAME, 'w')
    FF.write('\\\\ %s \n\n' % FnameTmp)
    if len(fba.objectives) > 0:
        if fba.objectives[fba.activeObjIdx].operation == None:
            print(
                '\nWARNING: Objective function \"{}\" has no \"operation\" defined assuming \"maximize\"'.format(
                    fba.objectives[fba.activeObjIdx].getId()
                )
            )
            fba.objectives[fba.activeObjIdx].operation = 'maximize'
            time.sleep(2)
        objO = fba.objectives[fba.activeObjIdx].operation.lower()
        objO = objO[0].upper() + objO[1:]
        FF.write('%s\n' % objO)
        objStr = '%s_objf: ' % fba.objectives[fba.activeObjIdx].getId()
        for fObj in fba.objectives[fba.activeObjIdx].flux_objectives:
            sign = None
            nc = 0.0
            try:
                nc = float(fObj.coefficient)
            except ValueError:
                if _HAVE_SYMPY_:
                    nc = sympy.Rational(fObj.coefficient).evalf()
                else:
                    raise ValueError(
                        'Invalid coefficient (%s) detected in fluxObjective %s'
                        % (fObj.coefficient, fObj.getId())
                    )
            if nc >= 0.0:
                sign = '+'
            else:
                sign = '-'
            if use_rational:
                objStr += ' %s %s%s%s' % (
                    sign,
                    sympy.Rational(abs(nc)),
                    multisymb,
                    fObj.reaction,
                )
            else:
                objStr += ' %s %s%s%s' % (sign, abs(nc), multisymb, fObj.reaction)
    else:
        objStr = '\n\\\\ No objectives defined\n\n'

    FF.write('%s\n' % objStr)
    CONST = BuildLPConstraints(fba, use_rational)
    CONST.seek(0)
    if fba.user_constraints != None and len(fba.user_constraints) > 0:
        UCONST = BuildLPUserConstraints(fba, use_rational=False)
        UCONST.seek(0)
    BOUNDS = BuildLPFluxBounds(fba, use_rational)
    BOUNDS.seek(0)
    if __DEBUG__:
        print(CONST.read())
        CONST.seek(0)
    if __DEBUG__:
        print(BOUNDS.read())
        BOUNDS.seek(0)

    FF.write('\nSubject To\n')
    FF.write(CONST.read())
    if fba.user_constraints != None and len(fba.user_constraints) > 0:
        FF.write('\\\\UserConstraints\n')
        FF.write(UCONST.read())
    FF.write('\nBounds\n')
    FF.write(BOUNDS.read())
    FF.write('\nEND\n')
    FF.close()
    if not quiet:
        print('writeModelLP has written a file to {}'.format(FNAME))
    return FNAME


def BuildHformatFluxBounds(fba, infinity_replace=None, use_rational=False):
    """
    Generate upper and lower flux bounds matrices for FBA optimization in H-format.

    This function constructs two matrices representing the upper and lower bounds
    of fluxes for reactions in a flux balance analysis (FBA) model. These matrices
    are suitable for linear programming (LP) formulations where flux bounds are
    applied as constraints in H-format.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model whose flux bounds will be represented in H-format matrices.
    infinity_replace : float or int, optional
        A value to replace ▒infinity in the flux bounds. If not provided, ▒infinity
        will be kept as is, which may not be suitable for all LP solvers.
    use_rational : bool, optional
        Whether to generate the matrices with rational numbers instead of floating-point
        numbers. This is useful for exact arithmetic in LP solvers that support rationals.

    Returns
    -------
    tuple of numpy.ndarray
        A tuple containing two matrices. The first matrix represents the lower bounds,
        and the second matrix represents the upper bounds for each reaction in the FBA
        model. The rows correspond to reactions, and the columns correspond to bounds.

    Example
    -------
    >>> from cbmpy import CBRead
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> lb_matrix, ub_matrix = BuildHformatFluxBounds(fba_model, infinity_replace=1000, use_rational=True)
    >>> print(lb_matrix)
    >>> print(ub_matrix)

    Notes
    -----
    H-format represents linear inequalities as Ax + b >= 0. The generated matrices need
    to be combined with the stoichiometry matrix and other constraints for a complete
    LP formulation of FBA.
    """
    LBs = {}
    UBs = {}
    for c in fba.flux_bounds:
        minv = None
        maxv = None
        oper = None
        R = c.reaction
        if c.operation in ['less', 'lessEqual']:
            if infinity_replace != None and numpy.isposinf([c.value])[0]:
                maxv = infinity_replace
            elif infinity_replace != None and numpy.isneginf([c.value])[0]:
                maxv = -infinity_replace
            else:
                maxv = c.value
            oper = c.operation
            if oper == 'less':
                oper = '<='
            else:
                oper = '<='
        if c.operation in ['greater', 'greaterEqual']:
            if infinity_replace != None and numpy.isposinf([c.value])[0]:
                minv = infinity_replace
            elif infinity_replace != None and numpy.isneginf([c.value])[0]:
                minv = -infinity_replace
            else:
                minv = c.value
            oper = c.operation
            if oper == 'greater':
                oper = '<='
            else:
                oper = '<='
        if c.operation == 'equals':
            if infinity_replace != None and numpy.isposinf([c.value])[0]:
                maxv = infinity_replace
            elif infinity_replace != None and numpy.isneginf([c.value])[0]:
                maxv = -infinity_replace
            else:
                maxv = c.value
            oper = '='

        if maxv != None:
            if not use_rational:
                UBs.update({R: float(maxv)})
            else:
                UBs.update({R: maxv})

        elif minv != None:
            if not use_rational:
                LBs.update({R: float(minv)})
            else:
                LBs.update({R: minv})

    BsRHS = []
    if not use_rational:
        LBm = numpy.zeros((len(LBs), fba.N.shape[1]))
        UBm = numpy.zeros((len(UBs), fba.N.shape[1]))
    else:
        if _HAVE_SYMPY_:
            LBm = sympy.zeros(len(LBs), fba.N.shape[1])
            UBm = sympy.zeros(len(UBs), fba.N.shape[1])
            print('\nUsing rational bounds')
        else:
            raise RuntimeError('\nError: SymPy required for rational operations')

    LBskeys = list(LBs)
    for lb in range(len(LBskeys)):
        if not use_rational:
            LBm[lb, fba.N.col.index(LBskeys[lb])] = 1.0
            BsRHS.append(LBs[LBskeys[lb]])
        else:
            LBm[lb, fba.N.col.index(LBskeys[lb])] = sympy.Rational(
                1.0
            ).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT'])
            BsRHS.append(
                sympy.Rational(LBs[LBskeys[lb]]).limit_denominator(
                    __CBCONFIG__['SYMPY_DENOM_LIMIT']
                )
            )

    UBskeys = list(UBs)
    for ub in range(len(UBskeys)):
        if not use_rational:
            UBm[ub, fba.N.col.index(UBskeys[ub])] = -1.0
            BsRHS.append(-UBs[UBskeys[ub]])
        else:
            UBm[ub, fba.N.col.index(UBskeys[ub])] = sympy.Rational(
                -1.0
            ).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT'])
            BsRHS.append(
                -sympy.Rational(UBs[UBskeys[ub]]).limit_denominator(
                    __CBCONFIG__['SYMPY_DENOM_LIMIT']
                )
            )

    del LBskeys, UBskeys

    return numpy.vstack([LBm, UBm]), BsRHS


def WriteModelHFormatFBA(
    fba,
    work_dir=None,
    use_rational=False,
    fullLP=True,
    format='%s',
    infinity_replace=None,
):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"\n'
    )
    time.sleep(1)
    writeModelHFormatFBA(fba, work_dir, use_rational, fullLP, format, infinity_replace)


def writeModelHFormatFBA(
    fba,
    work_dir=None,
    use_rational=False,
    fullLP=True,
    format='%s',
    infinity_replace=None,
):
    """
    Write an FBA-LP in polynomial H-Format file. This version has been replaced by `writeModelHFormatFBA2()`
    but is kept for backwards compatibility.

    Parameters
    ----------
    fba : CBMPy FBA model object
        The FBA model to write in H-format.
    work_dir : str, optional
        The directory where to save the output file. If not provided, uses the current directory.
    use_rational : bool, optional
        Whether to use rational numbers in the output. Requires sympy. Defaults to False.
    fullLP : bool, optional
        Whether to include the default objective function as a maximization target. Defaults to True.
    format : str, optional
        The format string for numbers. Defaults to '%s'.
    infinity_replace : float or None, optional
        A value to replace ▒infinity in bounds. If None, keeps ▒infinity as is.

    Returns
    -------
    str
        The filename of the written H-format file.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> fname = CBWrite.writeModelHFormatFBA(fba_model, work_dir='output/', use_rational=True)
    >>> print(fname)  # 'output/example_model_ine' or 'output/example_model_r.ine'

    """
    print('\nwriteModelHFormatFBA is deprecated please use: writeModelHFormatFBA2\n')

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')
    M = fba
    LHS = M.N.array.copy()
    RHS = [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
    LHS = numpy.vstack([LHS, -M.N.array.copy()])
    RHS += [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
    BsLHS, BsRHS = BuildHformatFluxBounds(M, infinity_replace=infinity_replace)
    if __DEBUG__:
        print(BsLHS)
        print(BsRHS)
    LHS = numpy.vstack([LHS, BsLHS])
    RHS += BsRHS
    del BsLHS, BsRHS
    if __DEBUG__:
        print(LHS)
        print(RHS)
    if not use_rational:
        name = M.getId().replace('.xml', '') + '.ine'
    else:
        name = M.getId().replace('.xml', '') + '_r.ine'

    RHS = numpy.array(RHS, 'd')
    RHS.shape = (len(RHS), 1)
    if __DEBUG__:
        print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    OBJ_FUNC = numpy.zeros(LHS.shape[1] + 1)
    for j in range(LHS.shape[1]):
        # first objective function, first flux objective
        if __DEBUG__:
            print(M.objectives[0].flux_objectives[0].reaction, M.N.col[j])
        if M.objectives[0].flux_objectives[0].reaction == M.N.col[j]:
            OBJ_FUNC[j] = float(M.objectives[0].flux_objectives[0].coefficient)
    if __DEBUG__:
        print(OBJ_FUNC)

    # for Ax >= B Hformat wants -B A >= 0
    LP = numpy.hstack([-RHS, LHS])
    OBJ_FUNC = numpy.hstack([-OBJ_FUNC[-1], OBJ_FUNC[:-1]])

    if __DEBUG__:
        print(OBJ_FUNC)
        print(LP)
    del LHS, RHS

    if work_dir == None:
        Fname = name
    else:
        assert os.path.exists(work_dir), '\nJanee ...'
        Fname = os.path.join(work_dir, name)

    F = open(Fname, 'w')
    F.write('* %s\n\n' % name)
    F.write('H-representation\n\nbegin\n')
    NUM_TYPE = 'real'
    if use_rational:
        NUM_TYPE = 'rational'
    F.write('%s  %s  %s\n' % (LP.shape[0], LP.shape[1], NUM_TYPE))

    strW = format + ' '
    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if not use_rational:
                if LP[r, c] == 0.0 or LP[r, c] == -0.0:
                    LP[r, c] = 0.0
                F.write(strW % LP[r, c])
            else:
                F.write('%s ' % sympy.Rational(format % LP[r, c]))
        F.write('\n')

    if fullLP:
        F.write('end\nlponly\n')
        F.write('maximize\n')
        for o in OBJ_FUNC:
            if not use_rational:
                F.write(strW % o)
            else:
                F.write('%s ' % sympy.Rational(format % o))
    else:
        F.write('end\n')
    F.write('\n')
    F.close()
    F = open(Fname.replace('.ine', '') + '.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return Fname


def WriteModelHFormatFBA2(
    fba,
    fname=None,
    work_dir=None,
    use_rational=False,
    fullLP=True,
    format='%s',
    infinity_replace=None,
):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"\n'
    )
    time.sleep(1)
    writeModelHFormatFBA2(
        fba, fname, work_dir, use_rational, fullLP, format, infinity_replace
    )


def writeModelHFormatFBA2(
    fba,
    fname=None,
    work_dir=None,
    use_rational=False,
    fullLP=True,
    format='%s',
    infinity_replace=None,
):
    """
    Write an FBA model to an H-Format file, suitable for use with linear programming solvers.

    This method is an enhanced version of `writeModelHFormatFBA()`, leveraging improvements
    such as the automatic use of rational numbers if a SymPy matrix is part of the input.

    Parameters
    ----------
    fba : PySCeS-CBM FBA object
        The FBA object to be written out.
    fname : str, optional
        The filename for the output file. If not specified the FBA object ID is used.
    work_dir : str, optional
        The directory where the file will be written. Defaults to the current working directory.
    use_rational : bool, optional
        If True, use rational numbers for output (requires Sympy). Defaults to False.
    fullLP : bool, optional
        If True, include the default objective function as a maximization target. Defaults to True.
    format : str, optional
        The format to use for numbers in the output. Defaults to '%s'.
    infinity_replace : float, optional
        A numeric value to replace any occurrences of infinity. If not specified, infinity values are left unchanged.

    Returns
    -------
    None

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> CBWrite.writeModelHFormatFBA2(fba_model, fname='model_h_format', use_rational=True)

    Note
    ----
    This method specifically prepares an H-Format file for further processing by LP solvers,
    allowing for the use of rational numbers for increased precision.
    """
    use_rational_old = use_rational
    if _HAVE_SYMPY_ and fba.N.__array_type__ == sympy.MutableDenseMatrix:
        use_rational = True
        print('INFO: using rational matrix')
    M = fba
    LHS = M.N.array.copy()
    # RHS = [0.0 for e in range(M.N.shape[0])]
    RHS = M.N.RHS.tolist()
    if __DEBUG__:
        print(LHS)
        print(RHS)
    LHS = numpy.vstack([LHS, -M.N.array.copy()])
    tmp = -M.N.RHS
    RHS += tmp.tolist()
    del tmp
    if __DEBUG__:
        print(LHS)
        print(RHS)
    BsLHS, BsRHS = BuildHformatFluxBounds(
        M, infinity_replace=infinity_replace, use_rational=use_rational
    )

    if use_rational:
        BsLHS = sympy.Matrix(BsLHS)
        den_lim = __CBCONFIG__['SYMPY_DENOM_LIMIT']
        BsRHS = [sympy.Rational(i).limit_denominator(den_lim) for i in BsRHS]
        for r in range(BsLHS.shape[0]):
            for c in range(BsLHS.shape[1]):
                BsLHS[r, c] = sympy.Rational(BsLHS[r, c]).limit_denominator(den_lim)

    if __DEBUG__:
        print(BsLHS)
        print(BsRHS)
    LHS = numpy.vstack([LHS, BsLHS])
    RHS += BsRHS
    del BsLHS, BsRHS
    if __DEBUG__:
        print(LHS)
        print(RHS)

    RHS = numpy.array(RHS)
    RHS.shape = (len(RHS), 1)
    if __DEBUG__:
        print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    OBJ_FUNC = numpy.zeros(LHS.shape[1] + 1)
    objIdx = M.activeObjIdx
    for j in range(LHS.shape[1]):
        for fo in range(len(M.objectives[objIdx].getFluxObjectiveReactions())):
            if M.objectives[objIdx].flux_objectives[fo].reaction == M.N.col[j]:
                print(M.objectives[objIdx].flux_objectives[fo].reaction, M.N.col[j])
                OBJ_FUNC[j] = float(M.objectives[objIdx].flux_objectives[fo].coefficient)
    ##  print OBJ_FUNC

    # for Ax >= B Hformat wants -B A >= 0
    LP = numpy.hstack([-RHS, LHS])
    OBJ_FUNC = numpy.hstack([-OBJ_FUNC[-1], OBJ_FUNC[:-1]])

    if __DEBUG__:
        print(OBJ_FUNC)
        print(LP)
    del LHS, RHS

    if work_dir != None:
        assert os.path.exists(work_dir), '\nJanee ...'
        fname = os.path.join(work_dir, fname)

    if fname == None:
        fname = M.getId().replace('.xml', '')
    if not use_rational:
        fname += '.ine'
    else:
        fname += '_r.ine'

    F = open(fname, 'w')
    F.write('* %s\n' % os.path.split(fname)[-1])
    F.write('H-representation\nbegin\n')
    NUM_TYPE = 'real'
    if use_rational:
        NUM_TYPE = 'rational'
    F.write('%s  %s  %s\n' % (LP.shape[0], LP.shape[1], NUM_TYPE))

    strW = format + ' '

    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if use_rational or use_rational_old:
                F.write('%s ' % LP[r, c])
            else:
                if LP[r, c] == 0.0 or LP[r, c] == -0.0:
                    LP[r, c] = 0.0
                F.write(strW % LP[r, c])
        F.write('\n')

    if fullLP:
        F.write('end\nlponly\n')
        F.write('maximize\n')  # check if Hformat has a minimize kw
        for o in OBJ_FUNC:
            if use_rational or use_rational_old:
                F.write('%s ' % sympy.Rational(format % o))
            else:
                if o == 0.0 or o == -0.0:
                    LP[r, c] = 0.0
                F.write(strW % o)
    else:
        F.write('end\n')
    F.write('\n')
    F.close()
    F = open(fname.replace('.ine', '') + '.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return fname


def writeStoichiometricMatrix(
    fba,
    fname=None,
    work_dir=None,
    use_rational=False,
    fullLP=True,
    format='%s',
    infinity_replace=None,
):
    """
    Write an FBA-LP in polynomial H-Format file. This is an improved version of `writeModelHFormatFBA()`
    which it replaces but is kept for backwards compatibility.

     - *fba*: PySCeS-CBM FBA object
        The FBA object to be written out.
     - *fname*: str, optional
        The output filename. If not defined, `fba.getId()` will be used.
     - *work_dir*: str, optional
        The output directory. If not defined, the current working directory is used.
     - *use_rational*: bool, optional
        Use rational numbers in output (requires SymPy). Defaults to False.
     - *fullLP*: bool, optional
        Include the default objective function as a maximization target. Defaults to True.
     - *format*: str, optional
        The number format string. Defaults to '%s'.
     - *infinity_replace*: numeric, optional
        A value to replace +/-inf in the bounds. If defined, this is the absolute value used.

    Returns
    -------
    str
        The filename where the stoichiometric matrix is saved.

    Example
    -------
    >>> fba_model = CBRead.readSBML3FBC('example_model.xml')
    >>> fname = writeStoichiometricMatrix(fba_model, use_rational=True)
    >>> print('Stoichiometric matrix saved to:', fname)

    """
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')
    M = fba
    LHS = M.N.array.copy()
    RHS = [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
    RHS += [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
    if __DEBUG__:
        print(BsLHS)
        print(BsRHS)
    if __DEBUG__:
        print(LHS)
        print(RHS)

    if __DEBUG__:
        print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    objIdx = M.activeObjIdx

    LP = LHS

    if __DEBUG__:
        print(OBJ_FUNC)
        print(LP)
    del LHS, RHS

    if work_dir != None:
        assert os.path.exists(work_dir), '\nJanee ...'
        fname = os.path.join(work_dir, fname)

    if fname == None:
        fname = M.getId().replace('.xml', '')
    if not use_rational:
        fname += '.ine'
    else:
        fname += '_r.ine'

    F = open(fname, 'w')

    strW = format + ' '
    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if not use_rational:
                if LP[r, c] == 0.0 or LP[r, c] == -0.0:
                    LP[r, c] = 0.0
                F.write(strW % LP[r, c])
            else:
                ##  print LP[r,c]
                F.write('%s ' % sympy.Rational(format % LP[r, c]))
        F.write('\n')

    F.close()
    F = open(fname.replace('.ine', '') + '.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return fname

def writeListToLP(
    fname, obj=None, const=None, bnds=None, work_dir=None, objtype='maximize'
):
    """
    Writes a linear programming model to a .lp file given objectives, constraints, and bounds.

    Parameters
    ----------
    fname : str
        The filename for the output file, without the extension.
    obj : list of str, optional
        A list of strings representing the objective function(s).
    const : list of str, optional
        A list of strings representing the constraints.
    bnds : list of str, optional
        A list of strings representing the bounds.
    work_dir : str, optional
        The directory where the .lp file will be saved. Defaults to the current working directory.
    objtype : {'maximize', 'minimize'}, default 'maximize'
        The optimization objective. Either 'maximize' or 'minimize'.

    Returns
    -------
    str
        The path to the created .lp file.

    Example
    -------
    >>> obj = ['obj: Z = 3x + y']
    >>> const = ['c1: -x + y <= 1', 'c2: x + y >= 2', 'c3: -2x + y <= 0']
    >>> bnds = ['x >= 0', 'y >= 0', 'x <= 4']
    >>> fname = 'example_model'
    >>> lp_file_path = writeListToLP(fname, obj=obj, const=const, bnds=bnds, objtype='maximize')
    >>> print(lp_file_path)

    """
    if work_dir is None:
        work_dir = os.getcwd()
    F = open(os.path.join(work_dir, fname + '.lp'), 'w')
    F.write("\\\\ %s\n" % fname)
    objtype = objtype.lower()
    if objtype == 'max':
        objtype = 'maximize'
    if objtype == 'min':
        objtype = 'minimize'
    if objtype in ['maximise', 'minimise']:
        objtype = objtype.replace('se', 'ze')
    assert objtype in ['maximize', 'minimize'], (
        "\nobjtype must be ['maximize', 'minimize'] not %s" % objtype
    )

    if obj is not None:
        F.write(f'\n{objtype.capitalize()}\n')
        for o in obj:
            F.write(f'{o}\n')
    if const is not None:
        F.write('\nSubject To\n')
        for c in const:
            F.write(f'{c}\n')
    if bnds is not None:
        F.write('\nBounds\n')
        for b in bnds:
            F.write(f'{b}\n')
    F.write('\nEnd\n')
    F.close()
    print(f'LP written to: {os.path.join(work_dir, fname)}.lp')
    return os.path.join(work_dir, fname + '.lp')


def writeMinDistanceLP(
    fname,
    fbas,
    work_dir=None,
    ignoreDistance=[],
    with_protein_cost=False,
    constraint_mode='strict',
    moma=False,
):
    """
    Write a linear programming (LP) problem to minimize the distance between multiple FBA models.

    Parameters
    ----------
    fname : str
        The filename for the output LP file.
    fbas : list of FBA model objects
        List of FBA model objects to be compared.
    work_dir : str, optional
        Directory where the LP file will be saved. If None, uses the current working directory.
    ignoreDistance : list of str, optional
        List of reaction IDs to be ignored in the distance calculation.
    with_protein_cost : bool, optional
        If True, includes protein cost in the calculation. Default is False.
    constraint_mode : str, optional
        Constraint mode to be used. Can be 'strict', 'math', or 'relaxed'. Default is 'strict'.
    moma : bool, optional
        If True, performs minimization of metabolic adjustment (MoMA) instead of FBA. Default is False.

    Returns
    -------
    str
        Path to the generated LP file.

    Example
    --------
    >>> fbas = [model1, model2]
    >>> lp_file = eMinDistanceLP('min_distance.lp', fbas, work_dir='./output', with_protein_cost=True)
    >>> print('LP file saved to:', lp_file)

    Notes
    -----
    The function generates an LP file that aims to minimize the flux distance
    between the provided FBA models, optionally considering the cost associated
    with protein expression. This functionality is particularly useful for comparing
    metabolic models or for integrating omics data.
    """
    if work_dir == None:
        work_dir = os.getcwd()

    fC = []
    objFname = ''
    for l in fbas:
        fC.append(len(l.reactions))
        objFname += l.prefix
    fC = numpy.array(fC)
    if __DEBUG__:
        print(fC)
        print((fC == fC[0]))
        print(numpy.alltrue((fC == fC[0])))
    if moma:
        assert len(fbas) == 2, '\nMOMA only defined for two inputs'
    else:
        assert numpy.alltrue(
            (fC == fC[0])
        ), '\nModels must have the same number of fluxes\n!'

    conL = []
    # model flux_bounds
    initial_cnstr = []

    if moma:
        initial_cnstr.append(BuildLPConstraintsStrict(fbas[0]))
        mt = fbas[1].clone()
        mt.setPrefix('mt_', 'all')
        initial_cnstr.append(BuildLPConstraintsStrict(mt))
        del mt
    else:
        for f in fbas:
            if constraint_mode == 'math':
                initial_cnstr.append(BuildLPConstraintsMath(f))
            elif constraint_mode == 'relaxed':
                initial_cnstr.append(BuildLPConstraintsRelaxed(f))
            else:
                initial_cnstr.append(BuildLPConstraintsStrict(f))

    for ib in initial_cnstr:
        ib.seek(0)
        for l in ib:
            conL.append(l.strip())
        conL.append(' ')
    del initial_cnstr

    bndL = []
    # model bounds
    initial_bnds = []
    for f in fbas:
        initial_bnds.append(BuildLPFluxBounds(f))
    for ib in initial_bnds:
        ib.seek(0)
        for l in ib:
            bndL.append(l.strip())
        bndL.append(' ')
    del initial_bnds

    artVar = []
    artVarX = []
    # ignoreDistance = []

    Combi = CBTools.ComboGen()
    Cnumber = 2
    Cdata = ''
    unique_combinations = None
    for x in range(len(fbas)):
        Cdata += '%s' % x
    if __DEBUG__:
        print(Cdata)

    Combi.uniqueCombinations(Cdata, Cnumber, temp=[])
    Combi.numberifyComb2Int()
    unique_combinations = Combi.combo_int

    if __DEBUG__:
        print('Data ({}):\n{}\n'.format(Cnumber, Cdata))
        print('UniqueCombinStr:\n{}'.format(Combi.combo))
        print('UniqueCombinations:\n{}'.format(unique_combinations))

    zbase = 0
    combcount = 1
    for uq in unique_combinations:
        if __DEBUG__:
            print(uq)
        MD0 = fbas[uq[0]]
        MD0 = fbas[uq[0]]
        if moma:
            MD1 = fbas[uq[1]].clone()
            MD1.setPrefix('mt_', 'all')
        else:
            MD1 = fbas[uq[1]]

        RiD1 = MD0.getReactionIds()
        RiD2 = MD1.getReactionIds()
        for s in range(len(MD0.reactions)):
            if __DEBUG__:
                print(RiD1[s], RiD2[s])
            if RiD1[s] not in ignoreDistance:
                if moma:
                    ridx2 = RiD2.index(MD1.prefix + RiD1[s])
                else:
                    ridx2 = s
                ##  av = 'z%s' % (zbase+s+1)
                av = 'zvar%s%s' % (combcount, RiD1[s].replace(MD0.prefix, ''))

                c1 = '%sa: %s - %s - %s <= 0.0' % (av, RiD1[s], RiD2[ridx2], av)
                ##  c1 = '%s - %s - %s <= 0.0' % (RiD1[s], RiD2[s], av)
                c2 = '%sb: %s - %s + %s >= 0.0' % (av, RiD1[s], RiD2[ridx2], av)
                ##  c2 = '%s - %s + %s >= 0.0' % (RiD1[s], RiD2[s], av)

                # add the protein cost
                if with_protein_cost:
                    av = '%s %s' % (MD0.reactions[s].annotation['CBM_PEPTIDE_COST'], av)
                artVar.append(av)
                conL.append(c1)
                conL.append(c2)

        combcount += 1
        zbase += len(MD0.reactions)

    objS = '%smulti: ' % objFname
    vcntr = 0
    for v in artVar:
        objS += '%s + ' % v
        vcntr += 1
        if vcntr >= 500:
            objS += '\n'
            vcntr = 0
    objS = objS[:-3]
    if len(artVarX) >= 1:
        objS += ' \\* Ignored: '
        for o in artVarX:
            objS += '%s ' % o
        objS += '*\\\n'
    objL = [objS]

    objFcnstr = [' ']
    assert (
        len(f.objectives[f.activeObjIdx].getFluxObjectiveReactions()) == 1
    ), "\nOnly single fluxObjectives dealt with at this time"
    if moma:
        f = fbas[0]
        objFcnstr.append(
            'C_%s: %s >= %f'
            % (
                f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                f.objectives[f.activeObjIdx].value,
            )
        )
    else:
        for f in fbas:
            objFcnstr.append(
                'C_%s: %s >= %f'
                % (
                    f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                    f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                    f.objectives[f.activeObjIdx].value,
                )
            )
    conL = conL + objFcnstr

    if __DEBUG__:
        print(objL)
        for c in conL:
            print(c)
        print(' ')
        for b in bndL:
            print(b)

    F = open(os.path.join(work_dir, fname + '.lp'), 'w')
    header = '\\\\ MultiInputMinimization: '
    for f in fbas:
        header += '%s, ' % f.getId()
    F.write('%s\n' % header[:-2])
    F.write('\nMinimize\n')
    for o in objL:
        F.write('%s \n' % o)
    F.write('\nSubject to\n')
    for c in conL:
        F.write('%s \n' % c)
    F.write('\nBounds\n')
    for b in bndL:
        F.write(' %s \n' % b)
    F.write('END\n\n')
    F.close()
    print('LP written to: {}.lp'.format(os.path.join(work_dir, fname)))
    return os.path.join(work_dir, fname + '.lp')


def writeMinDistanceLP_absL1(
    fname,
    fbas,
    work_dir=None,
    ignoreDistance=[],
    bigM=500,
    with_protein_cost=False,
    constraint_mode='strict',
    moma=False,
):
    """
    Writes a linear programming problem to minimize the absolute L1 distance between flux distributions of multiple FBA models.

    Parameters
    ----------
    fname : str
        Filename for the output .lp file.
    fbas : list of FBA model objects
        List of FBA models to compare.
    work_dir : str, optional
        Directory to write the output file. Defaults to current working directory.
    ignoreDistance : list of str, optional
        List of reaction IDs to ignore in the distance calculation.
    bigM : float, optional
        Big M constant used in linearization of the absolute value function.
    with_protein_cost : bool, optional
        If True, adds protein cost associated with flux to the optimization objective.
    constraint_mode : str, optional
        Indicates how constraints from FBA models are added ('strict', 'math', or 'relaxed').
    moma : bool, optional
        If True, sets MOMA (Minimization of Metabolic Adjusment) constraints instead of regular FBA constraints.

    Returns
    -------
    str
        Path to the generated LP file, relative to the working directory if specified.

    Example
    -------
    >>> fba_models = [model1, model2]
    >>> lp_file = writeMinDistanceLP_absL1("min_dist", fba_models, work_dir="outputs")
    >>> print("LP file written to:", lp_file)
    """
    if work_dir == None:
        work_dir = os.getcwd()

    fC = []
    objFname = ''
    for l in fbas:
        fC.append(len(l.reactions))
        objFname += l.prefix
    fC = numpy.array(fC)
    if __DEBUG__:
        print(fC)
        print((fC == fC[0]))
        print(numpy.alltrue((fC == fC[0])))
    if moma:
        assert len(fbas) == 2, '\nMOMA only defined for two inputs'
    else:
        assert numpy.alltrue(
            (fC == fC[0])
        ), '\nModels must have the same number of fluxes\n!'

    conL = []
    # model flux_bounds
    initial_cnstr = []

    if moma:
        initial_cnstr.append(BuildLPConstraintsStrict(fbas[0]))
        mt = fbas[1].clone()
        mt.setPrefix('mt_', 'all')
        initial_cnstr.append(BuildLPConstraintsStrict(mt))
        del mt
    else:
        for f in fbas:
            if constraint_mode == 'math':
                initial_cnstr.append(BuildLPConstraintsMath(f))
            elif constraint_mode == 'relaxed':
                initial_cnstr.append(BuildLPConstraintsRelaxed(f))
            else:
                initial_cnstr.append(BuildLPConstraintsStrict(f))

    for ib in initial_cnstr:
        ib.seek(0)
        for l in ib:
            conL.append(l.strip())
        conL.append(' ')
    del initial_cnstr

    bndL = []
    # model bounds
    initial_bnds = []
    bigMS = []
    for f in fbas:
        initial_bnds.append(BuildLPFluxBounds(f))
        ##  bigMS.append(max([abs(float(v.value)) for v in f.reactions]))
    for ib in initial_bnds:
        ib.seek(0)
        for l in ib:
            bndL.append(l.strip())
        bndL.append(' ')
    del initial_bnds

    Combi = CBTools.ComboGen()
    Cnumber = 2
    Cdata = ''
    unique_combinations = None
    for x in range(len(fbas)):
        Cdata += '%s' % x
    if __DEBUG__:
        print(Cdata)

    Combi.uniqueCombinations(Cdata, Cnumber, temp=[])
    Combi.numberifyComb2Int()
    unique_combinations = Combi.combo_int

    if __DEBUG__:
        print('Data ({}):\n{}\n'.format(Cnumber, Cdata))
        print('UniqueCombinStr:\n{}'.format(Combi.combo))
        print('UniqueCombinations:\n{}'.format(unique_combinations))
    print(bigMS)
    print('\nbigM = ', bigM, '\n')

    ##  ILPMETHOD = 'SK' # steven
    ILPMETHOD = 'GK'  # gunnar
    zbase = 0
    combcount = 1

    artVar = []
    artVarX = []
    # ignoreDistance = []
    boolVars = []

    for uq in unique_combinations:
        if __DEBUG__:
            print(uq)
        MD0 = fbas[uq[0]]
        if moma:
            MD1 = fbas[uq[1]].clone()
            MD1.setPrefix('mt_', 'all')
        else:
            MD1 = fbas[uq[1]]

        RiD1 = MD0.getReactionIds()
        RiD2 = MD1.getReactionIds()

        for s in range(len(MD0.reactions)):
            if __DEBUG__:
                print(RiD1[s], RiD2[s])
            if RiD1[s] not in ignoreDistance:
                if moma:
                    ridx2 = RiD2.index(MD1.prefix + RiD1[s])
                else:
                    ridx2 = s
                Var1 = RiD1[s]
                bVar1 = 'xvar_%s' % Var1
                absVar1 = 'absL_%s' % Var1
                Var2 = RiD2[ridx2]
                bVar2 = 'xvar_%s' % Var2
                absVar2 = 'absL_%s' % Var2

                c0a = '\n'
                c0a += '%s - %s >= 0\n' % (absVar1, Var1)
                c0a += '%s + %s >= 0\n' % (absVar1, Var1)

                if ILPMETHOD == 'GK':
                    # gunnar
                    c0a += '%s + %s - %s %s <= 0\n' % (absVar1, Var1, bigM, bVar1)
                    c0a += '%s - %s + %s %s <= %s\n' % (
                        absVar1,
                        Var1,
                        bigM,
                        bVar1,
                        bigM,
                    )
                elif ILPMETHOD == 'SK':
                    # steven
                    c0a += '%s - %s - %s %s <= 0\n' % (absVar1, Var1, bigM, bVar1)
                    c0a += '%s + %s + %s %s <= %s\n' % (
                        absVar1,
                        Var1,
                        bigM,
                        bVar1,
                        bigM,
                    )

                c0b = '\n'
                c0b += '%s - %s >= 0\n' % (absVar2, Var2)
                c0b += '%s + %s >= 0\n' % (absVar2, Var2)
                if ILPMETHOD == 'GK':
                    # gunnar
                    c0b += '%s + %s - %s %s <= 0\n' % (absVar2, Var2, bigM, bVar2)
                    c0b += '%s - %s + %s %s <= %s\n' % (
                        absVar2,
                        Var2,
                        bigM,
                        bVar2,
                        bigM,
                    )
                elif ILPMETHOD == 'SK':
                    # steven
                    c0b += '%s - %s - %s %s <= 0\n' % (absVar2, Var2, bigM, bVar2)
                    c0b += '%s + %s + %s %s <= %s\n' % (
                        absVar2,
                        Var2,
                        bigM,
                        bVar2,
                        bigM,
                    )

                c0 = c0a + c0b

                av = 'zvar%s%s' % (combcount, Var1.replace(MD0.prefix, ''))

                if bVar1 not in boolVars:
                    boolVars.append(bVar1)
                if bVar2 not in boolVars:
                    boolVars.append(bVar2)
                c1 = '%sa: %s - %s - %s <= 0.0' % (av, absVar1, absVar2, av)
                c2 = '%sb: %s - %s + %s >= 0.0' % (av, absVar1, absVar2, av)

                # add the protein cost
                if with_protein_cost:
                    av = '%s %s' % (MD0.reactions[s].annotation['CBM_PEPTIDE_COST'], av)
                artVar.append(av)

                conL.append(c0)
                conL.append(c1)
                conL.append(c2)
        combcount += 1
        zbase += len(MD0.reactions)

    objS = '%smulti: ' % objFname
    vcntr = 0
    for v in artVar:
        objS += '%s + ' % v
        vcntr += 1
        if vcntr >= 500:
            objS += '\n'
            vcntr = 0
    objS = objS[:-3]
    if len(artVarX) >= 1:
        objS += ' \\* Ignored: '
        for o in artVarX:
            objS += '%s ' % o
        objS += '*\\\n'
    objL = [objS]

    objFcnstr = [' ']
    assert (
        len(f.objectives[f.activeObjIdx].getFluxObjectiveReactions()) == 1
    ), "\nOnly single fluxObjectives dealt with at this time"
    if moma:
        # OFvalue = f.objectives[f.activeObjIdx].value
        f = fbas[0]
        objFcnstr.append(
            'C_%s: %s >= %f'
            % (
                f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                f.objectives[f.activeObjIdx].value,
            )
        )
    else:
        for f in fbas:
            objFcnstr.append(
                'C_%s: %s >= %f'
                % (
                    f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                    f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],
                    f.objectives[f.activeObjIdx].value,
                )
            )
    conL = conL + objFcnstr
    if __DEBUG__:
        print(objL)
        for c in conL:
            print(c)
        print(' ')
        for b in bndL:
            print(b)

    F = open(os.path.join(work_dir, fname + '.lp'), 'w')
    header = '\\\\ MultiInputMinimization: '
    for f in fbas:
        header += '%s, ' % f.getId()
    F.write('%s\n' % header[:-2])
    F.write('\nMinimize\n')
    for o in objL:
        F.write('%s \n' % o)
    F.write('\nSubject to\n')
    for c in conL:
        F.write('%s \n' % c)
    F.write('\nBounds\n')
    for b in bndL:
        F.write(' %s \n' % b)
    if len(boolVars) > 0:
        F.write('Binary\n')
        for b in boolVars:
            F.write(' %s \n' % b)
    F.write('\nEND\n')
    F.close()
    print('LP written to: {}.lp'.format(os.path.join(work_dir, fname)))
    return os.path.join(work_dir, fname + '.lp')


def writeMinDistanceLPwithCost(
    fname, fbas, work_dir=None, ignoreDistance=[], constraint_mode='strict'
):
    """
    For backwards compatability only
    """
    print(
        "\n\n**********\nDeprecation warning!\nPlease use writeMinDistanceLP(with_protein_cost=True) instead of writeMinDistanceLPwithCost()\n\n**********\n"
    )
    writeMinDistanceLP(
        fname,
        fbas,
        work_dir=work_dir,
        ignoreDistance=ignoreDistance,
        with_protein_cost=True,
        constraint_mode=constraint_mode,
    )


def writeOptimalSolution(fba, fname, Dir=None, separator=',', only_exchange=False):
    """
    Writes the optimal solution of a CBModel to a file, detailing reaction fluxes and gene associations.

    Parameters
    ----------
    fba : PySCeS CBModel object
        The CBModel instance containing the FBA optimization results.
    fname : str
        The base name for the output file. The file will be named `<fname>_solution.csv` or `<fname>_solution.txt` depending on the separator used.
    Dir : str, optional
        The directory where the output file will be saved. If not provided, uses the current directory.
    separator : str, optional
        The column separator in the output file. Defaults to ',', and the file extension will be `.csv`. If another separator is used, the file extension will be `.txt`.
    only_exchange : bool, optional
        If True, only reactions labeled as exchange reactions (those with `is_exchange` attribute set True) are included in the output.

    Returns
    -------
    None

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('path_to_model.xml')
    >>> CBWrite.writeOptimalSolution(model, 'model_optimal_solution', Dir='output_directory')

    """
    if Dir != None:
        assert os.path.exists(Dir), '\nPath does not exist'
        fname = os.path.join(Dir, fname)
    if separator == ',':
        fname_r = fname + '_solution.csv'
    else:
        fname_r = fname + '_solution.txt'
    objName = ''
    if len(fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions()) > 1:
        for J in fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions():
            objName += '%s_' % J
        objName = objName[:-1]
    else:
        objName = fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions()[0]
    try:
        F = open(fname_r, 'w')
    except IOError:
        print(
            '\nCSV file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(
                fname_r
            )
        )
        return
    cntr = 0
    F.write(
        '%s%s%s%s%s%s\n'
        % ('ObjectiveFunction', separator, objName, separator, separator, separator)
    )
    F.write(
        '\"%s\"%s%s%s%s%s%s%s%s%s\"%s\"%s\"%s\"\n'
        % (
            'Reaction',
            separator,
            'Value',
            separator,
            'LowerBound',
            separator,
            'UpperBound',
            separator,
            'Reduced cost',
            separator,
            'Name',
            separator,
            'Gene association',
        )
    )
    for r in fba.reactions:
        GO = False
        if not only_exchange:
            GO = True
        elif only_exchange and r.is_exchange:
            GO = True
        if GO:
            if 'GENE ASSOCIATION' in r.annotation:
                gene = r.annotation['GENE ASSOCIATION']
            else:
                gene = 'none'
            bnds = fba.getReactionBounds(r.getId())
            Lbnd = -numpy.inf
            Ubnd = numpy.inf
            if bnds != None:
                if bnds[1] != None:
                    Lbnd = bnds[1]
                if bnds[2] != None:
                    Ubnd = bnds[2]
            F.write(
                '\"%s\"%s%s%s%s%s%s%s%s%s\"%s\"%s\"%s\"\n'
                % (
                    r.getId(),
                    separator,
                    r.value,
                    separator,
                    Lbnd,
                    separator,
                    Ubnd,
                    separator,
                    r.reduced_cost,
                    separator,
                    r.name,
                    separator,
                    gene,
                )
            )
    F.flush()
    F.close()
    print('Reactions exported to {}'.format(fname_r))


def writeModelInfoToFile(
    fba, fname, Dir=None, separator=',', only_exchange=False, met_type='all'
):
    """
    Write comprehensive model information, including reactions and species, to file.

    Parameters
    ----------
    fba : CBModel
        The constraint-based model whose information will be written out.
    fname : str
        The filename for the output file.
    Dir : str, optional
        The directory where the output file will be saved. If not provided, the current working directory is used.
    separator : str, optional
        The column separator in the output file. Defaults to ','.
    only_exchange : bool, optional
        If True, only reactions labeled as exchange (those facilitating the exchange between the model and its environment) will be included in the output. Defaults to False.
    met_type : str, optional
        Specifies the type of species to include in the species output file. Can be 'all' for all species, 'boundary' for boundary species, or 'variable' for non-boundary species. Defaults to 'all'.

    Returns
    -------
    None

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('my_model.xml')
    >>> CBWrite.writeModelInfoToFile(model, 'model_info.csv')

    This will generate two files: 'model_info.rxns.csv' for reaction data and 'model_info.spec.csv' for species data.
    """
    writeReactionInfoToFile(
        fba, fname, Dir=Dir, separator=separator, only_exchange=only_exchange
    )
    writeSpeciesInfoToFile(fba, fname, Dir=Dir, separator=separator, met_type=met_type)


def writeReactionInfoToFile(fba, fname, Dir=None, separator=',', only_exchange=False):
    """
    Writes the reactions information from a CBModel to a file, optionally filtering only exchange reactions.

    Parameters
    ----------
    fba : CBModel
        The constraint-based model from which to write reaction information.
    fname : str
        The base name for the output file. The file will be named `<fname>.rxns.csv` if the separator is a comma,
        otherwise `<fname>.rxns.txt`.
    Dir : str, optional
        The directory where the output file will be saved. If not provided, the file is saved in the current directory.
    separator : str, optional
        The delimiter to use in the output file, defaults to ',' (comma), which results in a CSV file.
    only_exchange : bool, optional
        If True, only reactions marked as exchange reactions will be written to the file. Defaults to False.

    Raises
    ------
    IOError
        If the file cannot be opened for writing, possibly because it is locked by another application.

    Example
    -------
    >>> model = CBModel(...)
    >>> writeReactionInfoToFile(model, 'model_reactions', Dir='output/', only_exchange=True)
    >>> # This will write only exchange reactions to 'output/model_reactions.rxns.csv'.

    """
    if Dir != None:
        assert os.path.exists(Dir), '\nPath does not exist'
        fname = os.path.join(Dir, fname)
    if separator == ',':
        fname_r = fname + '.rxns.csv'
    else:
        fname_r = fname + '.rxns.txt'
    objName = ''
    if len(fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions()) > 1:
        for J in fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions():
            objName += '%s_' % J
        objName = objName[:-1]
    else:
        objName = fba.objectives[fba.activeObjIdx].getFluxObjectiveReactions()[0]
    try:
        F = open(fname_r, 'w')
    except IOError:
        print(
            '\nOutput file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(
                fname_r
            )
        )
        return
    cntr = 0
    F.write(
        '%s%s%s%s%s%s\n'
        % ('ObjectiveFunction', separator, objName, separator, separator, separator)
    )
    F.write(
        '\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n'
        % (
            'Reaction',
            separator,
            'LowerBound',
            separator,
            'UpperBound',
            separator,
            'Name',
            separator,
            'Equation',
            separator,
            'Gene association',
        )
    )
    for r in fba.reactions:
        if not only_exchange:
            GO = True
        elif only_exchange and r.is_exchange:
            GO = True
        if GO:
            if 'GENE ASSOCIATION' in r.annotation:
                gene = r.annotation['GENE ASSOCIATION']
            else:
                gene = 'none'
            bnds = fba.getReactionBounds(r.getId())
            Lbnd = -numpy.inf
            Ubnd = numpy.inf
            if bnds != None:
                if bnds[1] != None:
                    Lbnd = bnds[1]
                if bnds[2] != None:
                    Ubnd = bnds[2]
            if r.reversible:
                equation = ' %s ' % __CBCONFIG__['REVERSIBLE_SYMBOL']
            else:
                equation = ' %s ' % __CBCONFIG__['IRREVERSIBLE_SYMBOL']
            subs = ''
            prods = ''
            for rr in r.reagents:
                if rr.coefficient > 0.0:
                    if abs(rr.coefficient) == 1.0:
                        prods += ' + %s' % (rr.species_ref)
                    else:
                        prods += ' + %s %s' % (abs(rr.coefficient), rr.species_ref)
                else:
                    if abs(rr.coefficient) == 1.0:
                        subs += ' + %s' % (rr.species_ref)
                    else:
                        subs += ' + %s %s' % (abs(rr.coefficient), rr.species_ref)
            subs = subs[3:]
            prods = prods[3:]
            equation = subs + equation + prods
            F.write(
                '\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n'
                % (
                    r.getId(),
                    separator,
                    Lbnd,
                    separator,
                    Ubnd,
                    separator,
                    r.name,
                    separator,
                    equation,
                    separator,
                    gene,
                )
            )
            GO = False
    F.flush()
    F.close()
    print('Reactions exported to {}'.format(fname_r))


def writeSpeciesInfoToFile(fba, fname, Dir=None, separator=',', met_type='all'):
    """
    Writes comprehensive species information from a CBModel to file, including IDs, names, compartments, boundary status, chemical formulas, and charges.

    Parameters
    ----------
    fba : CBModel
        The constraint-based model whose species information will be written.
    fname : str
        The base filename for the output files. Two files will be generated, `<fname>.spec.<ext>` and `<fname>.s2r.<ext>`, where `<ext>` is either `csv` or `txt` based on `separator`.
    Dir : str, optional
        The directory where the output files will be saved. If provided and does not exist, an error is raised.
    separator : str, optional
        The separator to use in the output files (e.g., ',', ';', or '\t'). The default is ',' which generates CSV files.
    met_type : {'all', 'boundary', 'variable'}, optional
        Specifies which types of species to include in the output: 'all' includes all species, 'boundary' includes only boundary species, and 'variable' includes only non-boundary species. The default is 'all'.

    Raises
    ------
    IOError
        If the specified output files are locked or cannot be opened.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('path_to_model.xml')
    >>> CBWrite.writeSpeciesInfoToFile(model, 'model_species_info', Dir='./output', met_type='boundary')

    Notes
    -----
    This function generates two files per model: one detailing species information and another detailing the reactions in which each species is a reagent.
    """
    if Dir != None:
        assert os.path.exists(Dir), '\nPath does not exist'
        fname = os.path.join(Dir, fname)
    if separator == ',':
        fname_s = fname + '.spec.csv'
        fname_sr = fname + '.s2r.csv'
    else:
        fname_s = fname + '.spec.txt'
        fname_sr = fname + '.s2r.txt'
    if met_type not in ['all', 'boundary', 'variable']:
        met_type = 'all'

    try:
        F2 = open(fname_s, 'w')
        F3 = open(fname_sr, 'w')
    except IOError:
        print(
            '\nOutput file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(
                fname_s
            )
        )
        return

    maxReagentOf = 0
    F2.write(
        '\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n'
        % (
            'Species',
            separator,
            'Name',
            separator,
            'Compartment',
            separator,
            'Fixed',
            separator,
            'ChemicalFormula',
            separator,
            'Charge',
        )
    )
    for s in fba.species:
        if len(s.isReagentOf()) > maxReagentOf:
            maxReagentOf = len(s.reagent_of)
        GOS = False
        if met_type == 'all':
            GOS = True
        elif met_type == 'variable' and not s.is_boundary:
            GOS = True
        elif met_type == 'boundary' and s.is_boundary:
            GOS = True
        if GOS:
            F2.write(
                '\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n'
                % (
                    s.getId(),
                    separator,
                    s.getName(),
                    separator,
                    s.compartment,
                    separator,
                    s.is_boundary,
                    separator,
                    s.chemFormula,
                    separator,
                    s.charge,
                )
            )
    for s in fba.species:
        srl = len(s.reagent_of)
        rgOut = '\"%s\"%s' % (s.getId(), separator)

        GOS = False
        if met_type == 'all':
            GOS = True
        elif met_type == 'variable' and not s.is_boundary:
            GOS = True
        elif met_type == 'boundary' and s.is_boundary:
            GOS = True
        if GOS:
            for rg in range(maxReagentOf):
                if rg < srl:
                    rgOut += '%s%s' % (s.reagent_of[rg], separator)
                else:
                    rgOut += '%s' % separator
            F3.write('%s\n' % rgOut[:-1])
    F2.flush()
    F2.close()
    F3.flush()
    F3.close()
    print('Species exported to {}'.format(fname_s))
    print('Reagent map exported to {}'.format(fname_sr))


def printFBASolution(fba, include_all=False):
    """
    Prints the Flux Balance Analysis (FBA) optimal solution to the screen.

    This function prints the FBA optimal solution, including the objective function and its optimal value, to the screen.
    Optionally, it can also print all reaction fluxes in the FBA model.

    Parameters
    ----------
    fba : CBModel
        The CBModel instance containing the FBA optimization results.
    include_all : bool, optional
        If True, prints all reaction fluxes in the FBA model. Default is False, only printing the objective function's flux and value.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('path_to_model.xml')
    >>> CBWrite.printFBASolution(model, include_all=True)

    """
    OFflux = fba.objectives[fba.activeObjIdx].flux_objectives[0].reaction
    OFvalue = fba.objectives[fba.activeObjIdx].value
    OFSense = fba.objectives[fba.activeObjIdx].operation
    print('\n\n**********\nModel: {}\n\n'.format(fba.getId()))
    print('{} objective: {}\nOptimal value: {}\n\n'.format(OFSense, OFflux, OFvalue))
    if include_all:
        for J in fba.reactions:
            print('{}: {}'.format(J.getId(), J.value))
    print('**********\n')


def exportModel(fba, fname=None, fmt='lp', work_dir=None, use_rational='both'):
    """
    Export the FBA model in different formats.

    Parameters
    ----------
    fba : CBModel
        The FBA model object to export.
    fname : str, optional
        The filename for the exported model. If None, `fba.getId()` is used.
    fmt : str, optional
        The export format of the model. Options are 'lp' for CPLEX format, 'hformat' for Polyhedra format,
        or 'all' for both formats. Default is 'lp'.
    work_dir : str, optional
        The working directory where the exported file will be saved. If None, the current working directory is used.
    use_rational : str, optional
        Specifies if hformat files should be written using rational math. Applies if `fmt` is 'all' or 'hformat'.
        Can be 'both', True (rational math), or False (not using rational math). The default 'both' behavior is
        a legacy behavior and writes both types.

    Example
    --------
    >>> exportModel(fba_model, fname='fba_export', fmt='all', work_dir='./exports', use_rational='True')

    Notes
    -----
    The 'hformat' export ignores the 'fname' argument and only uses the model ID (`fba.getId()`) as a legacy behavior.

    """
    if work_dir == None:
        work_dir = os.getcwd()
    if fmt == 'all' or fmt == 'lp':
        writeModelLP(fba, work_dir=work_dir, fname=fname)
    if fmt == 'all' or fmt == 'hformat':
        if fname == None:
            fname = fba.getId().replace('.xml', '')
        if use_rational == 'both':
            writeModelHFormatFBA2(fba, fname=fname, work_dir=work_dir)
            writeModelHFormatFBA2(
                fba, fname=fname, work_dir=work_dir, use_rational=True
            )
        elif use_rational == True:
            writeModelHFormatFBA2(
                fba, fname=fname, work_dir=work_dir, use_rational=True
            )
        elif use_rational == False:
            writeModelHFormatFBA2(fba, fname=fname, work_dir=work_dir)


def writeProteinCostToCSV(fba, fname):
    """
    Writes the protein costs 'CBM_PEPTIDE_COST' annotation to a CSV file.

    Parameters
    ----------
    fba : CBModel
        An instantiated FBA object.
    fname : str
        The export file name.

    """
    F = open(fname + '.costs.csv', 'w')
    F.write('rid,minL,maxL,avgL,cost\n')
    for R in fba.reactions:
        rid = R.getId()
        pcost = ''
        avg_l = ''
        minL = ''
        maxL = ''
        if 'CBM_PEPTIDE_COST' in R.annotation:
            pcost = R.annotation['CBM_PEPTIDE_COST']
        if 'CBM_PEPTIDE_LENGTH_MAX' in R.annotation:
            maxL = R.annotation['CBM_PEPTIDE_LENGTH_MAX']
            if R.annotation['CBM_PEPTIDE_LENGTH_MAX'] == None:
                avg_l = 1
            else:
                if 'CBM_AVG_PEPTIDE_LENGTH' in R.annotation:
                    avg_l = R.annotation['CBM_AVG_PEPTIDE_LENGTH']
                else:
                    avg_l = 0
        if 'CBM_PEPTIDE_LENGTH_MIN' in R.annotation:
            minL = R.annotation['CBM_PEPTIDE_LENGTH_MIN']
        F.write('%s,%s,%s,%s,%s\n' % (rid, minL, maxL, avg_l, pcost))
    F.flush()
    F.close()
    print('Protein costs written to file: {}'.format(fname + '.csv'))


def WriteFVAtoCSV(id, fva, names, Dir=None, fbaObj=None):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV()\"
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV()\"\n'
    )
    time.sleep(1)
    writeFVAtoCSV(fva, names, id, Dir, fbaObj)


def writeFVAtoCSV(fvadata, names, fname, Dir=None, fbaObj=None):
    """
    Writes the result of a Flux Variability Analysis (FVA) method to a CSV file.

    This function takes as input the output from a Flux Variability Analysis and creates a
    CSV file detailing the optimal value, minimum and maximum fluxes, and additional information
    for each reaction analyzed.

    Parameters
    ----------
    fvadata : list of lists
        The OUTPUT_ARRAY from FluxVariabilityAnalysis(), each sublist contains FVA data for a reaction.
    names : list of str
        The OUTPUT_NAMES from FluxVariabilityAnalysis(), list of reaction names corresponding to `fvadata`.
    fname : str
        Base name for the exported CSV file. The file will be named `<fname>.fva.csv`.
    Dir : str, optional
        The directory where the CSV file will be saved. If not given, uses the current directory.
    fbaObj : CBModel object, optional
        If provided, extra model information (like reaction equation and gene association) is added to the output.

    """
    if Dir != None:
        Dir = os.path.join(Dir, fname + '.fva.csv')
    else:
        Dir = fname + '.fva.csv'
    F = open(Dir, 'w')
    if fbaObj == None:
        F.write('name,optval,min,max,diff,red cost,minstat,maxstat\n')
    else:
        F.write(
            'name,optval,min,max,diff,red cost,minstat,maxstat,"equation","subsystem","gene association","confidence level"\n'
        )
    for Jidx in range(len(names)):
        if names[Jidx] != None:
            name = names[Jidx]
            optval = fvadata[Jidx][0]
            rc = fvadata[Jidx][1]
            min = fvadata[Jidx][2]
            max = fvadata[Jidx][3]
            diff = fvadata[Jidx][4]
            minstat = fvadata[Jidx][5]
            maxstat = fvadata[Jidx][6]
            if fbaObj == None:
                F.write(
                    '%s,%s,%s,%s,%s,%s,%s,%s\n'
                    % (name, optval, min, max, diff, rc, minstat, maxstat)
                )
            else:
                xInf = []
                Ro = fbaObj.getReaction(name)
                for k in [
                    'Equation',
                    'SUBSYSTEM',
                    'GENE ASSOCIATION',
                    'Confidence Level',
                ]:
                    if k in Ro.annotation:
                        xInf.append(Ro.annotation[k])
                    else:
                        xInf.append('')
                F.write(
                    '%s,%s,%s,%s,%s,%s,%s,%s,"%s","%s","%s",%s\n'
                    % (
                        name,
                        optval,
                        min,
                        max,
                        diff,
                        rc,
                        minstat,
                        maxstat,
                        xInf[0],
                        xInf[1],
                        xInf[2],
                        xInf[3],
                    )
                )

    F.flush()
    F.close()
    print('FVA results written to: {}'.format(Dir))


def WriteFVAdata(
    fva,
    names,
    fname,
    work_dir=None,
    roundec=None,
    scale_min=False,
    appendfile=False,
    info=None,
):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeFVAdata()\"
    """
    print(
        '\nINFO: this method will be deprecated please update your scripts to use \"writeFVAdata()\"\n'
    )
    time.sleep(1)
    writeFVAdata(fva, names, fname, work_dir, roundec, scale_min, appendfile, info)


def writeFVAdata(
    fvadata,
    names,
    fname,
    work_dir=None,
    roundec=None,
    scale_min=False,
    appendfile=False,
    info=None,
):
    """
    Writes the results of a Flux Variability Analysis (FVA) to a CSV file.

    This function creates a CSV file containing columns for Flux Name, FVA Minimum, FVA Maximum,
    Optimal Value, and Span based on the results of an FVA analysis. It can optionally also
    append information to an existing file and include additional notes.

    Parameters
    ----------
    fvadata : list of lists
        The result of Flux Variability Analysis (FVA), typically OUTPUT_ARRAY from FluxVariabilityAnalysis().
    names : list of str
        The reaction names associated with `fvadata`, typically OUTPUT_NAMES from FluxVariabilityAnalysis().
    fname : str
        The base filename for the CSV output.
    work_dir : str, optional
        The output directory for the CSV files. If not provided, the current directory is used.
    roundec : int, optional
        The number of decimal places to round off the output. If `None`, no rounding is applied.
    scale_min : bool, default=False
        If `True`, normalize each flux such that the FVA_MIN value is 0.0.
    appendfile : bool, default=False
        If `True`, append the data to an existing file instead of creating a new one.
    info : str, optional
        A string added to the results as an extra column, useful with `appendfile`.

    Returns
    -------
    None

    Example
    -------
    >>> fvadata = [[0, 1, 2, 3], [0.1, -0.1, 0.5, 0.2]]
    >>> names = ['R1', 'R2']
    >>> writeFVAdata(fvadata, names, 'FVA_output', work_dir='results', roundec=2, appendfile=True, info='Additional Information')
    """
    if work_dir != None:
        work_dir = os.path.join(work_dir, fname + '.fvadata.csv')
    else:
        work_dir = fname + '.fvadata.csv'
    if not appendfile:
        F = open(work_dir, 'w')
    else:
        F = open(work_dir, 'a')
    if info == None:
        F.write('%s,%s,%s,%s,%s\n' % ('Jid', 'min', 'max', 'optval', 'span'))
    else:
        F.write('%s,%s,%s,%s,%s,%s\n' % ('Jid', 'min', 'max', 'optval', 'span', 'info'))
    for Jidx in range(len(names)):
        if names[Jidx] != None:
            name = names[Jidx]
            max = fvadata[Jidx][3]
            min = fvadata[Jidx][2]
            optval = fvadata[Jidx][0]
            if roundec != None:
                max = round(max, roundec)
                min = round(min, roundec)
                optval = round(optval, roundec)
            if scale_min:
                if min > 0.0:
                    max = max - min
                    optval = optval - min
                    min = 0.0
                elif min <= 0.0:
                    max = max + abs(min)
                    optval = optval + abs(min)
                    min = 0.0
            if info == None:
                if roundec == None:
                    F.write(
                        '%s,%s,%s,%s,%s\n'
                        % (
                            name,
                            min,
                            max,
                            optval,
                            abs(fvadata[Jidx][3] - fvadata[Jidx][2]),
                        )
                    )
                else:
                    F.write(
                        '%s,%s,%s,%s,%s\n'
                        % (
                            name,
                            min,
                            max,
                            optval,
                            round(abs(fvadata[Jidx][3] - fvadata[Jidx][2]), roundec),
                        )
                    )

            else:
                if roundec == None:
                    F.write(
                        '%s,%s,%s,%s,%s,%s\n'
                        % (
                            name,
                            min,
                            max,
                            optval,
                            abs(fvadata[Jidx][3] - fvadata[Jidx][2]),
                            info,
                        )
                    )
                else:
                    F.write(
                        '%s,%s,%s,%s,%s,%s\n'
                        % (
                            name,
                            min,
                            max,
                            optval,
                            round(abs(fvadata[Jidx][3] - fvadata[Jidx][2]), roundec),
                            info,
                        )
                    )
    F.flush()
    F.close()
    print('FVAdata results written to: {}'.format(work_dir))


def writeSolutions(fname, sols=[], sep=',', extra_output=None, fba=None):
    """
    Writes multiple flux balance analysis (FBA) solution flux dictionaries to a CSV file.

    This function exports the flux values of multiple FBA solutions to a CSV file. Each dictionary
    in the `sols` list should contain reaction IDs as keys and their corresponding flux values as values.

    Parameters
    ----------
    fname : str
        The path and file name where the solutions should be saved.
    sols : list of dict
        A list of dictionaries with reaction IDs as keys and flux values as values. Each dictionary represents a solution.
    sep : str, optional
        The column delimiter to use in the CSV file, by default ','.
    extra_output : CBModel, optional
        If provided, additional information from this CBModel is added to the output file, such as reaction names,
        by default None.
    fba : CBModel, optional
        The FBA model associated with the solutions, required if `extra_output` is True, by default None.

    Examples
    --------
    >>> sols = [{'R1': 10, 'R2': 5}, {'R1': -5, 'R2': 15}]
    >>> writeSolutions('output_solutions.csv', sols, sep=',', extra_output=None)

    The above example will create a file named 'output_solutions.csv' with the flux values from the `sols` list.
    If `extra_output` is given (a CBModel instance), reaction names and potentially other information
    will be included next to the reaction IDs.
    """
    assert len(sols) >= 2, "\nThere must be two or more solutions to work with"

    reac_ids = set([])
    for s in sols:
        reac_ids = reac_ids.union(set(s))
    reac_ids = list(reac_ids)
    reac_ids.sort()

    reac_names = []
    ##  reac_bnds = []
    ##  reac_eqns = []
    if extra_output != None and fba != None:
        for r in reac_ids:
            reac_names.append(fba.getReaction(r).getName())
            ##  reac_bnds.append(fba.getReactionBounds(r.getId()))

    F = open(fname + '.csv', 'w')
    for r in range(len(reac_ids)):
        row = '%s%s' % (reac_ids[r], sep)
        s_str = ''
        for s in range(len(sols)):
            if reac_ids[r] in sols[s]:
                s_str += '%s%s' % (sols[s][reac_ids[r]], sep)
            else:
                s_str += '%s%s' % ('\"none\"', sep)
        row += s_str
        if extra_output:
            row += '\"%s\"\n' % reac_names[r]
        else:
            row = row[:-1] + '\n'
        F.write(row)
    F.flush()
    F.close()
    print('\nSolutions written to: \"%s\"\n' % fname)


def generateBGID(num, prefix):
    """
    Create a BGID (BiGG ID) generator.

    This function creates a generator that yields a unique BiGG
    ID by incrementing a numerical part appended to a given prefix. The numerical
    part is zero-padded to ensure a fixed length.

    Parameters
    ----------
    num : int
        The starting number for the ID generation.
    prefix : str
        A two-letter prefix for the ID. Must be exactly two alphabetic characters.

    Yields
    ------
    str
        A string in the format `<prefix><num>` where `<num>` is zero-padded to 6 digits.

    Raises
    ------
    AssertionError
        If the provided prefix is not two alphabetic characters.

    Example
    -------
    >>> bgid_generator = generateBGID(0, 'ST')
    >>> next(bgid_generator)
    'ST000001'
    >>> next(bgid_generator)
    'ST000002'

    """
    assert prefix[0].isalpha() and prefix[1].isalpha(), '\nInvalid prefix'
    while True:
        num += 1
        yield '{}{:0>6}'.format(prefix, num)


def convertFloatToExcel(num, roundoff):
    """
    Converts a float to Excel compatible "number".

    This function handles special float values (NaN, positive infinity, negative infinity) by converting them
    into strings recognizable by Excel. Other numerical values are rounded to the specified number of digits and returned.

    Parameters
    ----------
    num : float
        The float number to be converted.
    roundoff : int
        The number of decimal places to round the number to.

    Returns
    -------
    str or float
        Returns a string for special float values ('', 'nan', 'inf', '-inf') or the rounded float number.

    Examples
    --------
    >>> convertFloatToExcel(0.123456789, 2)
    0.12

    >>> convertFloatToExcel(numpy.inf, 6)
    'inf'

    >>> convertFloatToExcel(numpy.nan, 3)
    'nan'

    >>> convertFloatToExcel(None, 0)
    ''

    Notes
    -----
    This function is kept in sync with `convertExcelToFloat` to ensure consistency between conversion operations.

    """
    ## keep this in sync with convertExcelToFloat
    if num is None:
        rval = ''
    elif numpy.isnan(num):
        rval = 'nan'
    elif numpy.isposinf(num):
        rval = 'inf'
    elif numpy.isneginf(num):
        rval = '-inf'
    else:
        rval = round(num, roundoff)
    return rval


def convertExcelToFloat(num):
    """
    Converts an Excel "number" to a float.

    Parameters
    ----------
    num : str
        The string representation of the number to convert. Special strings
        for float('nan'), float('+inf'), and float('-inf') are accepted.

    Returns
    -------
    float or None
        The converted float value, or None if input is an empty string.

    Examples
    --------
    >>> convertExcelToFloat('1.23')
    1.23
    >>> convertExcelToFloat('inf')
    inf
    >>> convertExcelToFloat('-inf')
    -inf
    >>> convertExcelToFloat('nan')
    nan
    >>> convertExcelToFloat('')
    None

    """
    ## keep this in sync with convertFloatToExcel
    if num == '':
        rval = None
    elif num == 'nan':
        rval = numpy.nan
    elif num == 'inf':
        rval = numpy.inf
    elif num == '-inf':
        rval = -numpy.inf
    else:
        rval = float(num)
    return rval


def writeModelToExcel97(fba, filename, roundoff=6):
    """
    Exports the model as an Excel 97 spreadsheet.

    This function exports the contents of a CBMPy model instance (`fba`)
    into an Excel 97 format spreadsheet. This includes various sheets detailing
    information, solution, reactions, metabolites, reaction network (by reactions and metabolites),
    MIRIAM annotations, compartments, and groups.

    Parameters
    ----------
    fba : CBModel
        The CBMPy model instance to export.
    filename : str
        The filename of the workbook to be created.
    roundoff : int, optional
        The number of decimal places to which numerical values are rounded. Default is 6.

    Raises
    ------
    IOError
        If the output file is locked by an external application or cannot be opened for writing.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite
    >>> model = CBRead.readSBML3FBC('example_model.xml')
    >>> CBWrite.writeModelToExcel97(model, 'model_export', roundoff=4)

    Notes
    -----
    Requires the `xlwt` package for writing Excel files. If `xlwt` is not available,
    an error message is printed and the function exits without writing the file.
    """

    if not _HAVE_XLWT_:
        print(
            '\nERROR: Cannot create Excel file, XLWT package not available (http://pypi.python.org/pypi/xlwt)'
        )
        return
    try:
        F = open('{}.xls'.format(filename), 'wb')
        F.close()
    except:
        print(
            '\nERROR: cannot open file "{}"! Please close workbook before writing!\n'.format(
                '{}.xls'.format(filename)
            )
        )
        return

    fontB = xlwt.Font()
    fontB.bold = True
    fontI = xlwt.Font()
    fontI.italic = True
    fontBI = xlwt.Font()
    fontBI.bold = True
    fontBI.italic = True
    styleBold = xlwt.XFStyle()
    styleBold.font = fontB
    styleBoldC = xlwt.XFStyle()
    styleBoldC.font = fontB
    styleBoldC.alignment.horz = xlwt.Alignment.HORZ_CENTER
    styleItalic = xlwt.XFStyle()
    styleItalic.font = fontI
    styleItalicC = xlwt.XFStyle()
    styleItalicC.font = fontI
    styleItalicC.alignment.horz = xlwt.Alignment.HORZ_CENTER
    styleBoldItalic = xlwt.XFStyle()
    styleBoldItalic.font = fontBI
    styleBoldItalicC = xlwt.XFStyle()
    styleBoldItalicC.font = fontBI
    styleBoldItalicC.alignment.horz = xlwt.Alignment.HORZ_CENTER
    styleHyper = xlwt.easyxf('font: underline single')
    styleHyper.alignment.horz = xlwt.Alignment.HORZ_CENTER

    wb = xlwt.Workbook(encoding='utf-8')

    wsInf = wb.add_sheet('info')
    wsSol = wb.add_sheet('solution')
    wsRe = wb.add_sheet('reactions')
    wsMet = wb.add_sheet('metabolites')
    wsSt = wb.add_sheet('network_react')
    wsStR = wb.add_sheet('network_metab')
    wsMiriam = wb.add_sheet('miriam')
    wsComp = wb.add_sheet('compartments')
    wsGrp = wb.add_sheet('groups')

    ridx = 0
    wsInf.write(ridx, 0, 'id', styleBold)
    wsInf.write(ridx, 1, fba.getId())
    ridx += 1
    wsInf.write(ridx, 0, 'name', styleBold)
    wsInf.write(ridx, 1, fba.getName())
    ridx += 1
    wsInf.write(ridx, 0, 'objective', styleBold)
    try:
        wsInf.write(ridx, 1, fba.getActiveObjective().getOperation())
    except:
        wsInf.write(ridx, 1, 'undefined')
    cidx = 2
    try:
        for t_ in fba.getActiveObjective().getFluxObjectiveData():
            try:
                wsInf.write(ridx, cidx, t_[0])
            except:
                wsInf.write(ridx, cidx, 'undefined')
            try:
                wsInf.write(ridx, cidx + 1, t_[1])
            except:
                wsInf.write(ridx, cidx + 1, 'undefined')
            cidx += 2
    except:
        pass
    ridx += 1
    wsInf.write(ridx, 0, 'obj. value', styleBold)
    try:
        wsInf.write(ridx, 1, fba.getActiveObjective().getValue())
    except:
        wsInf.write(ridx, 1, 'undefined')
    ridx += 1
    wsInf.write(ridx, 0, 'modelcreated', styleBold)
    if fba.DATE_CREATED != None:
        for t_ in range(len(fba.DATE_CREATED)):
            wsInf.write(ridx, t_ + 1, fba.DATE_CREATED[t_])
    ridx += 1
    wsInf.write(ridx, 0, 'lastmodified', styleBold)
    if fba.DATE_MODIFIED != None:
        for t_ in range(len(fba.DATE_MODIFIED)):
            wsInf.write(ridx, t_ + 1, fba.DATE_MODIFIED[t_])
    ridx += 1
    cTime = list(time.gmtime())
    wsInf.write(ridx, 0, 'exported', styleBold)
    for t_ in range(len(cTime)):
        if t_ <= 5:
            wsInf.write(ridx, t_ + 1, cTime[t_])
    ridx += 1
    ccridx = 0
    for cr_ in fba.getModelCreators():
        cidx = 0
        ccridx += 1
        wsInf.write(ridx, cidx, 'creator{}'.format(ccridx), styleBold)
        for ne_ in cr_:
            cidx += 1
            wsInf.write(ridx, cidx, ne_)
        ridx += 1
    ridx += 1
    wsInf.write(
        ridx,
        0,
        xlwt.Formula(
            'HYPERLINK("{}";"{}")'.format(
                'http://cbmpy.sourceforge.net',
                'Exported by CBMPy {}'.format(__version__),
            )
        ),
        styleHyper,
    )
    ridx += 1


def writeModelToCOMBINEarchive(
    mod,
    fname=None,
    directory=None,
    sbmlname=None,
    withExcel=True,
    vc_given='CBMPy',
    vc_family='Software',
    vc_email='None',
    vc_org='cbmpy.sourceforge.net',
    add_cbmpy_annot=True,
    add_cobra_annot=True,
):
    """
    Writes a model in SBML and Excel format to a COMBINE archive.

    This function exports a CBMPy model to a COMBINE archive containing the model in SBML format
    and optionally in a human-readable Excel format. Additionally, it includes metadata
    describing the model and the archive creation details.

    Parameters
    ----------
    mod : CBModel
        The CBMPy model instance to export.
    fname : str, optional
        The base filename of the output archive. The final archive will be `<fname>.sbex.zip`.
        If not given, must be specified before calling the function.
    directory : str, optional
        The directory where the output archive will be saved. If not specified, the current
        directory is used.
    sbmlname : str, optional
        The filename for the SBML export inside the archive. If not specified, defaults
        to `<fname>.xml`.
    withExcel : bool, optional
        If True (default), includes an Excel version of the model in the archive.
    vc_given : str, optional
        The given name (e.g., first name) of the creator for metadata purposes. Defaults to 'CBMPy'.
    vc_family : str, optional
        The family name (e.g., surname) of the creator for metadata purposes. Defaults to 'Software'.
    vc_email : str, optional
        The email of the creator for metadata purposes. Defaults to 'None'.
    vc_org : str, optional
        The organization of the creator for metadata purposes. Defaults to 'cbmpy.sourceforge.net'.
    add_cbmpy_annot : bool, optional
        If True (default), adds CBMPy KeyValueData annotations to the SBML export.
    add_cobra_annot : bool, optional
        If True (default), adds COBRA-style <notes> annotations to the SBML export.

    Raises
    ------
    AssertionError
        If the temporary directory for assembling the archive could not be created.

    Example
    -------
    >>> from cbmpy import CBRead, CBWrite, CBModel
    >>> model = CBRead.readSBML3FBC('path_to_your_model_file.xml')
    >>> CBWrite.writeModelToCOMBINEarchive(model, fname='model_archive', directory='output_path',
    ...                                    sbmlname='model_sbml', withExcel=True)

    Notes
    -----
    The generated COMBINE archive (.sbex.zip) is compatible with tools that understand the COMBINE
    archive specification and includes metadata in RDF format according to the OMEX Metadata guidelines.

    """
    scTime = time.strftime('%Y-%m-%dT%H:%M:%S') + '%i:00' % (time.timezone / 60 / 60)
    if directory != None:
        zfpath = os.path.join(directory, fname + '.sbex.zip')
    else:
        zfpath = fname + '.sbex.zip'
    zf = zipfile.ZipFile(zfpath, mode='w', compression=zipfile.ZIP_DEFLATED)
    if sbmlname != None:
        if sbmlname.endswith('.xml'):
            sbmlf = sbmlname
            xlf = sbmlname[:-4]
        else:
            sbmlf = sbmlname + '.xml'
            xlf = sbmlname
    else:
        sbmlf = fname + '.xml'
        xlf = fname
    ptmp = os.path.join(os.getcwd(), 'sedxtmp')
    if not os.path.exists(ptmp):
        os.makedirs(ptmp)
    assert os.path.exists(
        ptmp
    ), "Could not create temporary archive directory: {}".format(ptmp)

    MFstr = ''
    MDstr = ''
    MFstr += '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">\n'
    MFstr += ' <content location="." format="http://identifiers.org/combine.specifications/omex"/>\n'
    MFstr += ' <content location="./metadata.rdf" format="http://identifiers.org/combine.specifications/omex-metadata"/>\n'

    # SBML
    writeSBML3FBCV2(
        mod,
        sbmlf,
        ptmp,
        add_cbmpy_annot=add_cbmpy_annot,
        add_cobra_annot=add_cobra_annot,
    )
    zf.write(os.path.join(ptmp, sbmlf), arcname=sbmlf)
    MFstr += ' <content location="./{}" format="http://identifiers.org/combine.specifications/sbml.level-3.version-1"/>\n'.format(
        sbmlf
    )

    # Excel
    if withExcel and _HAVE_XLWT_:
        writeModelToExcel97(mod, os.path.join(ptmp, xlf))
        xlf += '.xls'
        zf.write(os.path.join(ptmp, xlf), arcname=xlf)
        MFstr += ' <content location="./{}" format="http://mediatypes.appspot.com/application/vnd.ms-excel"/>'.format(
            xlf
        )
    MF = open(os.path.join(ptmp, 'manifest.xml'), 'w')
    MF.write(
        '<?xml version="1.0" encoding="utf-8"?>\n{}\n</omexManifest>\n'.format(MFstr)
    )
    MF.close()

    MD = open(os.path.join(ptmp, 'metadata.rdf'), 'w')
    MD.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    MD.write('<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
    MD.write('    xmlns:dcterms="http://purl.org/dc/terms/"\n')
    MD.write('    xmlns:vCard="http://www.w3.org/2006/vcard/ns#"\n')
    MD.write('    xmlns:bqmodel="http://biomodels.net/models-qualifiers">\n')
    MD.write(' <rdf:Description rdf:about=".">\n')
    # MDstr += '   <dcterms:description>\n     %s\n    </dcterms:description>\n' % self.omex_description
    MDstr += ' <dcterms:creator>\n'
    MDstr += ' <rdf:Bag>\n'
    MDstr += '  <rdf:li rdf:parseType="Resource">\n'
    MDstr += '   <vCard:hasName rdf:parseType="Resource">\n'
    MDstr += '    <vCard:family-name>{}</vCard:family-name>\n'.format(vc_family)
    MDstr += '    <vCard:given-name>{}</vCard:given-name>\n'.format(vc_given)
    MDstr += '   </vCard:hasName>\n'
    MDstr += '   <vCard:hasEmail rdf:resource="{}" />\n'.format(vc_email)
    MDstr += '   <vCard:organization-name>\n'
    MDstr += '      {}\n'.format(vc_org)
    MDstr += '   </vCard:organization-name>\n'
    MDstr += '  </rdf:li>\n'
    MDstr += ' </rdf:Bag>\n'
    MDstr += ' </dcterms:creator>\n'
    MDstr += '   <dcterms:created rdf:parseType="Resource">\n'
    MDstr += '    <dcterms:W3CDTF>{}</dcterms:W3CDTF>\n'.format(scTime)
    MDstr += '   </dcterms:created>\n'
    MDstr += '   <dcterms:modified rdf:parseType="Resource">\n'
    MDstr += '    <dcterms:W3CDTF>{}</dcterms:W3CDTF>\n'.format(scTime)
    MDstr += '   </dcterms:modified>\n'
    MD.write('{}'.format(MDstr))
    MD.write(' </rdf:Description>\n')
    MD.write('</rdf:RDF> \n')
    MD.close()

    zf.write(os.path.join(ptmp, 'manifest.xml'), arcname='manifest.xml')
    zf.write(os.path.join(ptmp, 'metadata.rdf'), arcname='metadata.rdf')
    zf.close()

    for f_ in os.listdir(ptmp):
        os.remove(os.path.join(ptmp, f_))
    try:
        os.removedirs(ptmp)
    except WindowsError:
        pass
    print('COMBINE archive created: {}'.format(fname + '.zip'))
