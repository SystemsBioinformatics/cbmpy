"""
CBMPy: CBWrite module
=====================
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
Last edit: $Author: bgoli $ ($Id: CBWrite.py 648 2018-05-23 20:01:54Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

if 'cDir' in vars():
    cDir=vars()['cDir']
else:
    cDir=None
import os, time, numpy, zipfile
# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio

from . import CBTools, CBXML, CBDataStruct


from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

_HAVE_SYMPY_ = False
try:
    import sympy
    if int(sympy.__version__.split('.')[1]) >= 7 and int(sympy.__version__.split('.')[2]) >= 5:
        _HAVE_SYMPY_ = True
    elif int(sympy.__version__.split('.')[0]) >= 1:
        _HAVE_SYMPY_ = True
    else:
        del sympy
        print('\nWARNING: SymPy version 0.7.5 or newer is required for symbolic matrix support.')
except ImportError:
    _HAVE_SYMPY_ = False
    print('\nERROR: SymPy import error (required for symbolic matrix support only).')

_HAVE_XLWT_ = False
try:
    import xlwt
    _HAVE_XLWT_ = True
except ImportError:
    print('\nINFO: No xlwt module available, Excel spreadsheet creation disabled')


def saveModel(model, filename, compress=False):
    """
    Saves the model to an SBML file using the lates SBML3 FBC version.
    
    - *model* the CBMPy model
    - *filename* the filename to write
    
    """
    
    writeSBML3FBCV2(model, filename, directory=None, gpr_from_annot=False, add_groups=True,\
                    add_cbmpy_annot=True, add_cobra_annot=False, validate=False, compress_bounds=True,\
                    zip_model=compress, return_model_string=False)
    

def writeSBML3FBC(fba, fname, directory=None, gpr_from_annot=False,\
                       add_groups=True, add_cbmpy_annot=True, add_cobra_annot=False,\
                       xoptions={'fbc_version': 1, 'validate' : False, 'compress_bounds' : True}):
    """
    Takes an FBA model object and writes it to file as SBML L3 FBC:

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname
     - *gpr_from_annot* [default=True] if enabled will attempt to add the gene protein associations from the annotations
       if no gene protein association objects exist
     - *add_groups* [default=True] add SBML3 groups (if supported by libSBML)
     - *add_cbmpy_annot* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
     - *add_cobra_annot* [default=True] add COBRA <notes> annotation
     - *xoptions* extended options

       - *fbc_version* [default=1] write SBML3FBC using version 1 (2013) or version 2 (2015)
       - *validate* [default=False] validate the output SBML file
       - *compress_bounds* [default=False] try compress output flux bound parameters
       - *zip_model* [default=False] compress the model using PKZIP encoding
       - *return_model_string* [default=False] return the SBML XML file as a string



    """
    sbml_level_version = (3,1)
    autofix=True
    return_fbc=False
    return CBXML.sbml_writeSBML3FBC(fba, fname, directory, sbml_level_version, autofix, return_fbc,\
                                    gpr_from_annot, add_groups, add_cbmpy_annot, add_cobra_annot, xoptions)


def writeSBML3FBCV2(fba, fname, directory=None, gpr_from_annot=False, add_groups=True, add_cbmpy_annot=True, add_cobra_annot=False,\
                    validate=False, compress_bounds=True, zip_model=False, return_model_string=False):
    """
    Takes an FBA model object and writes it to file as SBML L3 FBCv2 :

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname
     - *gpr_from_annot* [default=False] if enabled will attempt to add the gene protein associations from the annotations
     - *add_groups* [default=True] add SBML3 groups (if supported by libSBML)
     - *add_cbmpy_annot* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
     - *add_cobra_annot* [default=False] add COBRA <notes> annotation
     - *validate* [default=False] validate the output SBML file
     - *compress_bounds* [default=True] try compress output flux bound parameters
     - *zip_model* [default=False] compress the model using ZIP encoding
     - *return_model_string* [default=False] return the SBML XML file as a string

    """

    xoptions = {'fbc_version': 2, 'validate' : validate, 'compress_bounds' : compress_bounds, 'return_model_string' : return_model_string, 'zip_model' : zip_model}
    sbml_level_version=(3,1)
    autofix=True,
    return_fbc=False
    #if fbc_version == 2:
        #add_cobra_annot = False
    return CBXML.sbml_writeSBML3FBC(fba, fname, directory, sbml_level_version, autofix, return_fbc,\
                                    gpr_from_annot, add_groups, add_cbmpy_annot, add_cobra_annot, xoptions)


def writeCOBRASBML(fba, fname, directory=None):
    """
    Takes an FBA model object and writes it to file as a COBRA compatible :

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname

    """
    return CBXML.sbml_writeCOBRASBML(fba, fname, directory)

def writeSBML2FBA(fba, fname, directory=None, sbml_level_version=None):
    """
    Takes an FBA model object and writes it to file as SBML L2 with FBA annotations.
    Note if you want to write BiGG/FAME style annotations then you must use *sbml_level_version=(2,1)*

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *sbml_level_version* [default=None] a tuple containing the SBML level and version e.g. (2,1)

    This is a utility wrapper for the function `CBXML.sbml_writeSBML2FBA`

    """
    CBXML.sbml_writeSBML2FBA(fba, fname, directory, sbml_level_version)


def writeSensitivitiesToCSV(sensitivities, fname):
    """
    Write out a sensitivity report using the objective sensitivities and
    bound sensitivity dictionaries created by e.g. cplx_getSensitivities().

     - *sensitivity* tuple containing

      - *obj_sens* dictionary of objective coefficient sensitivities (per flux)
      - *rhs_sens* dictionary of constraint rhs sensitivities (per constraint)
      - *bound_sens* dictionary of bound sensitivities (per flux)

     - *fname* output filename e.g. fname.csv

    """
    obj_sens = sensitivities[0]
    rhs_sens = sensitivities[1]
    bound_sens = sensitivities[2]
    F = open(fname+'_flux_sensitivity.csv', 'w')
    head = "Flux,Reduced cost,OCS low,OC value,OCS high,LB low,LB high,UB low,UB high"
    F.write(head+'\n')
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
        F.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (j,rc,lcs,ocv,ucs,lbs,lb,ub,ubs))

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
            F.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (j,rc,lcs,ocv,ucs,lbs,lb,ub,ubs))
    F.flush()
    F.close()

    F = open(fname+'_constraint_sensitivity.csv', 'w')
    F.write('Constraint,RHS low,RHS,RHS high\n')
    for c in rhs_sens:
        F.write('%s,%s,%s,%s\n' % (c, rhs_sens[c][0], rhs_sens[c][1], rhs_sens[c][2]))
    F.flush()
    F.close()

def WriteModelRaw(fba, work_dir=cDir):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelRaw()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeModelRaw()\"\n')
    time.sleep(1)
    writeModelRaw(fba, work_dir)

def writeModelRaw(fba, work_dir=cDir):
    """
    Writes a fba (actually just dumps it) to a text file.

     - *fba* an instantiated FBAmodel instance
     - *work_dir* directory designated for output

    """
    if work_dir == None:
        work_dir = os.getcwd()
    FF = open(os.path.join(work_dir,'WriteModelRawOutput.txt'), 'w')
    FF.write('Species information:\n\n')
    for s in fba.species:
        FF.write('%s: value=%f, is_boundary=%s, name=%s\n' % (s.getId(), s.value, int(s.is_boundary), s.name))
    FF.write('\nReaction information:\n\n')
    for r in fba.reactions:
        FF.write('%s: reversible=%s, name=%s\n' % (r.getId(), int(r.reversible), r.name))
        FF.write('\t%s\n' % r.getSpeciesIds())
    FF.write('\nConstraint information:\n\n')
    for c in fba.flux_bounds:
        FF.write('%s %s %f\n' % (c.reaction, c.operation, c.value))
    FF.write('\nObjective information:\n\n')
    FF.write('Active Objective: %s (%s)\n' % (fba.objectives[fba.activeObjIdx].getId(), fba.objectives[fba.activeObjIdx].operation))
    for o in fba.objectives:
        FF.write('%s: %s, %s\n' % (o.getId(), o.operation, o.getFluxObjectiveReactions()))
    if hasattr(fba.N, 'shape'):
        FF.write('\nStoichiometric information:\n\n')
        FF.write('N-matrix dimensions = (%s,%s)\n' % fba.N.shape)
    FF.close()
    print('WriteModelRaw has written a file to {}'.format(os.path.join(work_dir,'WriteModelRawOutput.txt')))

def BuildLPFluxBounds(fba, use_rational=False):
    """
    Build and return a csio that contains the flux bounds in LP format
    """
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    FFS = csio.StringIO()
    c2s = {}
    for c in fba.flux_bounds:
        if __DEBUG__: print('%s: %s' % (c.getId(), c.operation))
        minv = None
        maxv = None
        oper = None
        R = c.reaction
        if c.operation in ['less','lessEqual']:
            maxv = c.value
            oper = c.operation
            if oper == 'less':
                oper = '<='
                if __DEBUG__: print('LP Bounds operator must be <= not <')
            else:
                oper = '<='
        if c.operation in ['greater','greaterEqual']:
            minv = c.value
            oper = c.operation
            if oper == 'greater':
                oper = '<='
                if __DEBUG__: print('LP Bounds operator must be <= not <')
            else:
                oper = '<='
        if c.operation == 'equals':
            maxv = c.value
            oper = '='
        if maxv != None:
            if use_rational:
                c2s.update({R+'r' : '%s %s' % (oper, sympy.Rational(maxv))})
            else:
                if maxv == float('inf'):
                    c2s.update({R+'r' : '%s +%s' % (oper, maxv)})
                else:
                    c2s.update({R+'r' : '%s %s' % (oper, maxv)})
        elif minv != None:
            if use_rational:
                c2s.update({R+'l' : '%s %s' % (sympy.Rational(minv), oper)})
            else:
                c2s.update({R+'l' : '%s %s' % (minv, oper)})

    c2sk = list(c2s)
    c2sk.sort()
    for r in (c2sk):
        rht = None
        lht = None
        R = r[:-1]
        if R+'r' in c2s:
            rht = c2s[R+'r']
            c2s.pop(R+'r')
        if R+'l' in c2s:
            lht = c2s[R+'l']
            c2s.pop(R+'l')
        if lht != None and rht != None:
            if __DEBUG__: print('{} {} {}'.format(lht, R, rht))
            ##  FFS.write('%s: %s %s %s\n' % (R, lht, R, rht))
            FFS.write('%s %s %s\n' % (lht, R, rht)) # THIS MUST BE SO for GLPK
        elif lht == None and rht == None:
            if __DEBUG__: print('Skipping: ({}, {}, {})'.format(lht, R, rht))
        elif lht != None:
            if __DEBUG__: print('{} {}'.format(lht, R))
            #FFS.write('%s: %s %s\n' % (R, lht, R))
            FFS.write('%s %s\n' % (lht, R))
        elif rht != None:
            if __DEBUG__: print('{} {}'.format(R, rht))
            #FFS.write('%s: %s %s\n' % (R, R, rht))
            FFS.write('%s %s\n' % (R, rht))
        else:
            print('CONFUSION: ({}, {}, {})'.format(lht, R, rht))
    return FFS

def BuildLPConstraints(fba, use_rational=False):
    """
    Build and return a csio that contains constraint constructed from
    the StoichiometeryLP object

     - *fba* an fba model object which has a stoichiometry
     - *use_rational* write rational number output [default=False]
    """

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational output')
    rebuild_stoich = False
    if not hasattr(fba,'N') or fba.N == None:
        rebuild_stoich = True
    else:
        if len([s for s in fba.species if not s.is_boundary]) != fba.N.array.shape[0]:
            rebuild_stoich = True
        elif len(fba.reactions) != fba.N.array.shape[1]:
            rebuild_stoich = True
    if rebuild_stoich:
        print('\nWarning FBA object has inconsistant stoichiometric matrix, rebuilding it now.')
        CBTools.addStoichToFBAModel(fba)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r,c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName : RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    ##  for r in constrsk:
    for r in range(len(fba.N.row)):
        if len(constr[fba.N.row[r]]) > 0:
            FFS.write(' %s: ' % fba.N.row[r])
            ##  FFS.write(' ')
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
    Build and return a csio that contains constraint constructed from
    the StoichiometeryLP object

     - *fba* an fba model object which has a stoichiometry
     - *use_rational* write rational number output [default=False]
    """

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational output')
    rebuild_stoich = False
    assert fba.user_constraints != None, "\nNo user constraints to build"
    if not hasattr(fba,'CM') or fba.CM == None:
        rebuild_stoich = True

    if rebuild_stoich:
        print('\nWarning FBA object has inconsistant user constraint matrix, rebuilding it now.')
        CBTools.addStoichToFBAModel(fba)

    constr = {}
    for r in range(fba.CM.array.shape[0]):
        rowName = fba.CM.row[r]
        RCon = []
        for c in range(fba.CM.array.shape[1]):
            colName = fba.CM.col[c]
            colCoef = fba.CM.array[r,c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName : RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    ##  for r in constrsk:
    for r in range(len(fba.CM.row)):
        if len(constr[fba.CM.row[r]]) > 0:
            FFS.write(' %s: ' % fba.CM.row[r])
            ##  FFS.write(' ')
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
    Build and return a csio that contains the constaints in LP format
    Relaxed refers to dS/dt >= 0
    """
    raise DeprecationWarning("\nThis method is deprecated")
    if not hasattr(fba,'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r,c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName : RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write('%s: ' % r)
            ##  FFS.write(' ')
            for col in constr[r]:
                if col[0] > 0.0:
                    FFS.write('+%.20f %s ' % (col[0], col[1]))
                else:
                    FFS.write('%.20f %s ' % (col[0], col[1]))
            FFS.write('>= 0\n')

    return FFS

def BuildLPConstraintsStrict(fba, use_rational=False):
    """
    Build and return a csio that contains the constaints in LP format
    Strict refers to dS/dt = 0
    """
    ##  print "Consider using the new BuildLPConstraints() method"

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    if not hasattr(fba,'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r,c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName : RCon})

    FFS = csio.StringIO()
    constrsk = list(constr)
    constrsk.sort()
    ##  for r in constrsk:
    for r in fba.N.row:
        if len(constr[r]) > 0:
            FFS.write(' %s: ' % r)
            ##  FFS.write(' ')
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
    Build and return a csio that contains the constaints in LP format
    Strict refers to dS/dt => 0 and dS/dt <= 0
    """
    raise DeprecationWarning("\nThis method is deprecated")
    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('Warning: install Sympy for rational IO')

    if not hasattr(fba,'N') or fba.N == None:
        print('\nWarning FBA object has no stoichiometric matrix constructing it now.')
        CBTools.addStoichToFBAModel(fba)
        time.sleep(1)

    constr = {}
    for r in range(fba.N.array.shape[0]):
        rowName = fba.N.row[r]
        RCon = []
        for c in range(fba.N.array.shape[1]):
            colName = fba.N.col[c]
            colCoef = fba.N.array[r,c]
            if colCoef != 0.0:
                RCon.append((colCoef, colName))
        constr.update({rowName : RCon})

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

def WriteModelLPOld(fba, work_dir=None, multisymb=' ', lpt=True, constraint_mode='strict', use_rational=False, format='%s'):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelLPOld()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeModelLPOld()\"\n')
    time.sleep(1)
    writeModelLPOld(fba, work_dir, multisymb, lpt, constraint_mode, use_rational, format)

def writeModelLPOld(fba, work_dir=None, multisymb=' ', lpt=True, constraint_mode='strict', use_rational=False, format='%s'):
    """
    Writes a fba as an LP/LPT

     - *fba* an instantiated FBAmodel instance
     - *work_dir* directory designated for output
     - *multisymb* the multiplication symbol (default: <space>)
     - *lpt* the file format (default: True for lpt) or False for lp

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
        FnameTmp = FnameTmp+'.rat'
    if not lpt:
        FNAME = FnameTmp+'.lp'
        FF = open(FNAME, 'w')
        FF.write('Problem\n %s\n\n' % FnameTmp)
    else:
        FNAME = FnameTmp+'.lp'
        FF = open(FNAME, 'w')
        FF.write('\\\\ %s \n\n' % FnameTmp)
    objO = fba.objectives[fba.activeObjIdx].operation.lower()
    objO = objO[0].upper() + objO[1:]
    FF.write('%s\n' % objO)
    objStr = '%s_objf: ' % fba.objectives[fba.activeObjIdx].getId()
    for fObj in fba.objectives[fba.activeObjIdx].fluxObjectives:
        sign = None
        nc = 0.0
        try:
            nc = float(fObj.coefficient)
        except ValueError:
            print('Suspected rational number ({}) detected in fluxObjective {}'.format(fObj.coefficient, fObj.getId()))
        if nc >= 0.0:
            sign = '+'
        else:
            sign = '-'
        # TODO: if use_rational is not used simply try and evaluate the coefficient string with
        # sympy.Rational.evalf()  and use this as the value for nc ... also remove use_rational case
        if use_rational:
            objStr += ' %s%s%s' % (sympy.Rational(fObj.coefficient), multisymb, fObj.reaction)
            ##  FF.write('%s: %s%s%s\n' % (fobj0.reaction, sympy.Rational(fobj0.coefficient), multisymb, fobj0.reaction))
        else:
            objStr += ' %s %s%s%s' % (sign, abs(nc), multisymb, fObj.reaction)
            ##  FF.write('%s: %s%s%s\n' % (fobj0.reaction, fobj0.coefficient, multisymb, fobj0.reaction))
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
    if __DEBUG__: print(CONST.read()); CONST.seek(0)
    if __DEBUG__: print(BOUNDS.read()); BOUNDS.seek(0)

    FF.write('\nSubject To\n')
    FF.write(CONST.read())
    FF.write('\nBounds\n')
    FF.write(BOUNDS.read())
    FF.write('\nEND\n')
    FF.close()
    print('writeModelLP has written a file to {}'.format(NAME))
    return FNAME

def WriteModelLP(fba, work_dir=None, fname=None, multisymb=' ', format='%s', use_rational=False, constraint_mode=None, quiet=False):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelLP()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeModelLP()\"\n')
    time.sleep(1)
    writeModelLP(fba, work_dir, fname, multisymb, format, use_rational, constraint_mode, quiet)

def writeModelLP(fba, work_dir=None, fname=None, multisymb=' ', format='%s', use_rational=False, constraint_mode=None, quiet=False):
    """
    Writes an FBA object as an LP in CPLEX LP format

     - *fba* an instantiated FBAmodel instance
     - *work_dir* directory designated for output
     - *fname* the file name [default=fba.getId()]
     - *multisymb* the multiplication symbol (default: <space>)
     - *format* the number format of the output
     - *use_rational* output rational numbers [default=False]
     - *quiet* [default=False] supress information messages

    """
    if constraint_mode != None:
        print("\nConstraint_mode has been deprecated")
        time.sleep(5)

    if not _HAVE_SYMPY_ and use_rational:
        use_rational = False
        print('\nWarning switching to floating point arithmetic: install Sympy for rational IO')
        time.sleep(2)

    FNAME = None
    if fname == None:
        fname = fba.getId()
    if work_dir != None:
        FnameTmp = os.path.join(work_dir,  fname)
    else:
        FnameTmp = fname

    FNAME = FnameTmp+'.lp'
    FF = open(FNAME, 'w')
    FF.write('\\\\ %s \n\n' % FnameTmp)
    if len(fba.objectives) > 0:
        ##  print fba.objectives
        if fba.objectives[fba.activeObjIdx].operation == None:
            print('\nWARNING: Objective function \"{}\" has no \"operation\" defined assuming \"maximize\"'.format(fba.objectives[fba.activeObjIdx].getId()))
            fba.objectives[fba.activeObjIdx].operation = 'maximize'
            time.sleep(2)
        objO = fba.objectives[fba.activeObjIdx].operation.lower()
        objO = objO[0].upper() + objO[1:]
        FF.write('%s\n' % objO)
        objStr = '%s_objf: ' % fba.objectives[fba.activeObjIdx].getId()
        for fObj in fba.objectives[fba.activeObjIdx].fluxObjectives:
            sign = None
            nc = 0.0
            try:
                nc = float(fObj.coefficient)
            except ValueError:
                if _HAVE_SYMPY_:
                    nc = sympy.Rational(fObj.coefficient).evalf()
                else:
                    raise ValuError( 'Invalid coefficient (%s) detected in fluxObjective %s' % (fObj.coefficient, fObj.getId()))
            if nc >= 0.0:
                sign = '+'
            else:
                sign = '-'
            if use_rational:
                objStr += ' %s %s%s%s' % (sign, sympy.Rational(abs(nc)), multisymb, fObj.reaction)
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
    if __DEBUG__: print(CONST.read()); CONST.seek(0)
    if __DEBUG__: print(BOUNDS.read()); BOUNDS.seek(0)

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
    Build and return a csio that contains the flux bounds in H format

     - *fba* a PySCeS-CBM FBA object
     - *infinity_replace* [default=None] if defined this is the abs(value) of +-<infinity>

    """
    LBs = {}
    UBs = {}
    for c in fba.flux_bounds:
        ##  print '%s: %s' % (c.getId(), c.operation)
        minv = None
        maxv = None
        oper = None
        R = c.reaction
        if c.operation in ['less','lessEqual']:
            if infinity_replace != None and numpy.isposinf([c.value])[0]:
                maxv = infinity_replace
            elif infinity_replace != None and numpy.isneginf([c.value])[0]:
                maxv = -infinity_replace
            else:
                maxv = c.value
            ##  print 'maxv', maxv
            oper = c.operation
            if oper == 'less':
                oper = '<='
                ##  print 'LP Bounds operator must be <= not <'
            else:
                oper = '<='
        if c.operation in ['greater','greaterEqual']:
            if infinity_replace != None and numpy.isposinf([c.value])[0]:
                minv = infinity_replace
            elif infinity_replace != None and numpy.isneginf([c.value])[0]:
                minv = -infinity_replace
            else:
                minv = c.value
            ##  print 'minv', minv
            oper = c.operation
            if oper == 'greater':
                oper = '<='
                ##  print 'LP Bounds operator must be <= not <'
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
            ##  print 'maxv2', maxv
        if maxv != None:
            if not use_rational:
                UBs.update({R : float(maxv)})
            else:
                UBs.update({R : maxv})

        elif minv != None:
            if not use_rational:
                LBs.update({R : float(minv)})
            else:
                LBs.update({R : minv})

    if __DEBUG__:
        print(' ')
        print(LBs)
        print(UBs)

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
            LBm[lb, fba.N.col.index(LBskeys[lb])] = sympy.Rational(1.0).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT'])
            print(LBs[LBskeys[lb]])
            BsRHS.append(sympy.Rational(LBs[LBskeys[lb]]).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT']))

    if __DEBUG__:
        print(fba.N.col)
        print(LBm)
        print(UBm)
        print(BsRHS)

    UBskeys = list(UBs)
    for ub in range(len(UBskeys)):
        if not use_rational:
            UBm[ub, fba.N.col.index(UBskeys[ub])] = -1.0
            BsRHS.append(-UBs[UBskeys[ub]])
        else:
            UBm[ub, fba.N.col.index(UBskeys[ub])] = sympy.Rational(-1.0).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT'])
            BsRHS.append(-sympy.Rational(UBs[UBskeys[ub]]).limit_denominator(__CBCONFIG__['SYMPY_DENOM_LIMIT']))

    del LBskeys, UBskeys

    if __DEBUG__:
        print(fba.N.col)
        print(LBm)
        print(UBm)
        print(BsRHS)
        print(' ')

    return numpy.vstack([LBm, UBm]), BsRHS

def WriteModelHFormatFBA(fba, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"\n')
    time.sleep(1)
    writeModelHFormatFBA(fba, work_dir, use_rational, fullLP, format, infinity_replace)

def writeModelHFormatFBA(fba, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    Write an FBA-LP in polynomial H-Format file. This version has been replaced by `writeModelHFormatFBA2()`
    but is kept for backwards compatability.

     - *fba* a PySCeS-CBM FBA object
     - *Work_dir* [default=None] the output directory
     - *use_rational* [default=false] use rational numbers in output (requires sympy)
     - *fullLP* [default=True] include the default objective function as a maximization target
     - *format* [default='%s'] the number format string
     - *infinity_replace* [default=None] if defined this is the abs(value) of +-<infinity>

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

    RHS = numpy.array(RHS,'d')
    RHS.shape = (len(RHS), 1)
    if __DEBUG__: print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    OBJ_FUNC = numpy.zeros(LHS.shape[1]+1)
    for j in range(LHS.shape[1]):
        # first objective function, first flux objective
        if __DEBUG__: print(M.objectives[0].fluxObjectives[0].reaction, M.N.col[j])
        if M.objectives[0].fluxObjectives[0].reaction == M.N.col[j]:
            OBJ_FUNC[j] = float(M.objectives[0].fluxObjectives[0].coefficient)
    if __DEBUG__: print(OBJ_FUNC)

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

    strW = format+' '
    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if not use_rational:
                if LP[r,c] == 0.0 or LP[r,c] == -0.0:
                    LP[r,c] = 0.0
                F.write(strW % LP[r,c])
            else:
                F.write('%s ' % sympy.Rational(format % LP[r,c]))
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
    F = open(Fname.replace('.ine','')+'.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return Fname

def WriteModelHFormatFBA2(fba, fname=None, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeModelHFormatFBA2()\"\n')
    time.sleep(1)
    writeModelHFormatFBA2(fba, fname, work_dir, use_rational, fullLP, format, infinity_replace)

def writeModelHFormatFBA2(fba, fname=None, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    Write an FBA-LP in polynomial H-Format file. This is an improved version of `WriteModelHFormatFBA()`
    which it replaces. Note that if a SymPy matrix is used as input then use_rational is automatically enabled.

     - *fba* a PySCeS-CBM FBA object
     - *fname* [default=None] the output filename, fba.getId() if not defined
     - *Work_dir* [default=None] the output directory
     - *use_rational* [default=false] use rational numbers in output (requires sympy)
     - *fullLP* [default=True] include the default objective function as a maximization target
     - *format* [default='%s'] the number format string
     - *infinity_replace* [default=None] if defined this is the abs(value) of +-<infinity>

    """
    use_rational_old = use_rational
    if _HAVE_SYMPY_ and fba.N.__array_type__ == sympy.MutableDenseMatrix:
        use_rational = True
        print('INFO: using rational matrix')
    M = fba
    LHS = M.N.array.copy()
    #RHS = [0.0 for e in range(M.N.shape[0])]
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
    BsLHS, BsRHS = BuildHformatFluxBounds(M, infinity_replace=infinity_replace, use_rational=use_rational)

    if use_rational:
        BsLHS = sympy.Matrix(BsLHS)
        den_lim = __CBCONFIG__['SYMPY_DENOM_LIMIT']
        BsRHS = [sympy.Rational(i).limit_denominator(den_lim) for i in BsRHS]
        for r in range(BsLHS.shape[0]):
            for c in range(BsLHS.shape[1]):
                BsLHS[r,c] = sympy.Rational(BsLHS[r,c]).limit_denominator(den_lim)


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
    if __DEBUG__: print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    OBJ_FUNC = numpy.zeros(LHS.shape[1]+1)
    objIdx = M.activeObjIdx
    for j in range(LHS.shape[1]):
        for fo in range(len(M.objectives[objIdx].getFluxObjectiveReactions())):
            if M.objectives[objIdx].fluxObjectives[fo].reaction == M.N.col[j]:
                print(M.objectives[objIdx].fluxObjectives[fo].reaction, M.N.col[j])
                OBJ_FUNC[j] = float(M.objectives[objIdx].fluxObjectives[fo].coefficient)
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

    strW = format+' '

    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if use_rational or use_rational_old:
                F.write('%s ' % LP[r,c])
            else:
                if LP[r,c] == 0.0 or LP[r,c] == -0.0:
                    LP[r,c] = 0.0
                F.write(strW % LP[r,c])
        F.write('\n')

    if fullLP:
        F.write('end\nlponly\n')
        F.write('maximize\n') # check if Hformat has a minimize kw
        for o in OBJ_FUNC:
            #if not use_rational:
                #F.write(strW % o)
            #else:
                #F.write('%s ' % sympy.Rational(format % o))

            if use_rational or use_rational_old:
                F.write('%s ' % sympy.Rational(format % o))
            else:
                if o == 0.0 or o == -0.0:
                    LP[r,c] = 0.0
                F.write(strW % o)

        # then we can use use this
        ##  F.write('%s\n' % M.objectives[M.activeObjIdx].operation)
        ##  if M.activeObjIdx].operation == 'maximize':
            ##  for o in OBJ_FUNC:
                ##  if not use_rational:
                    ##  F.write(strW % o)
                ##  else:
                    ##  F.write('%s ' % sympy.Rational(format % o))
        ##  else:
            ##  for o in OBJ_FUNC:
                ##  o = -o
                ##  if not use_rational:
                    ##  F.write(strW % o)
                ##  else:
                    ##  F.write('%s ' % sympy.Rational(format % o))
    else:
        F.write('end\n')
    F.write('\n')
    F.close()
    F = open(fname.replace('.ine','')+'.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return fname


def writeStoichiometricMatrix(fba, fname=None, work_dir=None, use_rational=False, fullLP=True, format='%s', infinity_replace=None):
    """
    Write an FBA-LP in polynomial H-Format file. This is an improved version of `WriteModelHFormatFBA()`
    which it replaces but is kept for backwards compatability.

     - *fba* a PySCeS-CBM FBA object
     - *fname* [default=None] the output filename, fba.getId() if not defined
     - *Work_dir* [default=None] the output directory
     - *use_rational* [default=false] use rational numbers in output (requires sympy)
     - *fullLP* [default=True] include the default objective function as a maximization target
     - *format* [default='%s'] the number format string
     - *infinity_replace* [default=None] if defined this is the abs(value) of +-<infinity>

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
    #LHS = numpy.vstack([LHS, -M.N.array.copy()])
    RHS += [0.0 for e in range(M.N.shape[0])]
    if __DEBUG__:
        print(LHS)
        print(RHS)
    #BsLHS, BsRHS = BuildHformatFluxBounds(M, infinity_replace=infinity_replace)
    if __DEBUG__:
        print(BsLHS)
        print(BsRHS)
    #LHS = numpy.vstack([LHS, BsLHS])
    #RHS += BsRHS
    #del BsLHS, BsRHS
    if __DEBUG__:
        print(LHS)
        print(RHS)

    #RHS = numpy.array(RHS,'d')
    #RHS.shape = (len(RHS), 1)
    if __DEBUG__: print(RHS)
    ##  LP = numpy.hstack([LHS, RHS])

    #OBJ_FUNC = numpy.zeros(LHS.shape[1]+1)
    objIdx = M.activeObjIdx
    #for j in range(LHS.shape[1]):
        #for fo in range(len(M.objectives[objIdx].getFluxObjectiveReactions())):
            #if M.objectives[objIdx].fluxObjectives[fo].reaction == M.N.col[j]:
                #print(M.objectives[objIdx].fluxObjectives[fo].reaction, M.N.col[j])
                #OBJ_FUNC[j] = float(M.objectives[objIdx].fluxObjectives[fo].coefficient)
    ###  print OBJ_FUNC

    # for Ax >= B Hformat wants -B A >= 0
    #LP = numpy.hstack([-RHS, LHS])
    #OBJ_FUNC = numpy.hstack([-OBJ_FUNC[-1], OBJ_FUNC[:-1]])

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
    #F.write('* %s\n' % os.path.split(fname)[-1])
    #F.write('H-representation\nbegin\n')
    #NUM_TYPE = 'real'
    #if use_rational:
        #NUM_TYPE = 'rational'
    #F.write('%s  %s  %s\n' % (LP.shape[0], LP.shape[1], NUM_TYPE))

    strW = format+' '
    for r in range(LP.shape[0]):
        for c in range(LP.shape[1]):
            if not use_rational:
                if LP[r,c] == 0.0 or LP[r,c] == -0.0:
                    LP[r,c] = 0.0
                F.write(strW % LP[r,c])
            else:
                ##  print LP[r,c]
                F.write('%s ' % sympy.Rational(format % LP[r,c]))
        F.write('\n')

    #if fullLP:
        #F.write('end\nlponly\n')
        #F.write('maximize\n') # check if Hformat has a minimize kw
        #for o in OBJ_FUNC:
            #if not use_rational:
                #F.write(strW % o)
            #else:
                #F.write('%s ' % sympy.Rational(format % o))
        # then we can use use this
        ##  F.write('%s\n' % M.objectives[M.activeObjIdx].operation)
        ##  if M.activeObjIdx].operation == 'maximize':
            ##  for o in OBJ_FUNC:
                ##  if not use_rational:
                    ##  F.write(strW % o)
                ##  else:
                    ##  F.write('%s ' % sympy.Rational(format % o))
        ##  else:
            ##  for o in OBJ_FUNC:
                ##  o = -o
                ##  if not use_rational:
                    ##  F.write(strW % o)
                ##  else:
                    ##  F.write('%s ' % sympy.Rational(format % o))
    #else:
        #F.write('end\n')
    #F.write('\n')
    F.close()
    F = open(fname.replace('.ine','')+'.columns.txt', 'w')
    for j in range(M.N.array.shape[1]):
        F.write('%s,%s\n' % (j, M.N.col[j]))
    F.write('\n')
    F.close()
    return fname


def writeListToLP(fname, obj=None, const=None, bnds=None, work_dir=None, objtype='maximize'):
    if work_dir == None:
        work_dir = os.getcwd()
    F = open(os.path.join(work_dir, fname+'.lp'), 'w')
    F.write("\\\\ %s\n" % fname)
    objtype = objtype.lower()
    if objtype == 'max': objtype = 'maximize'
    if objtype == 'min': objtype = 'minimize'
    if objtype in ['maximise', 'minimise']:
        objtype = objtype.replace('se','ze')
    assert objtype in ['maximize', 'minimize'], "\nobjtype must be ['maximize', 'minimize'] not %s" % objtype

    if obj != None:
        if objtype == 'maximize':
            F.write('\nMaximize\n')
        elif objtype == 'minimize':
            F.write('\nMinimize\n')
        for o in obj:
            F.write(' %s\n' % o)
    if const != None:
        F.write('\nSubject to\n')
        for c in const:
            F.write('%s\n' % c)
    if bnds != None:
        F.write('\nBounds\n')
        for b in bnds:
            F.write('%s\n' % b)
    F.write('\nEND\n')
    F.close()
    print('LP written to: {}.lp'.format(os.path.join(work_dir, fname)))
    return os.path.join(work_dir, fname+'.lp')


def writeMinDistanceLP(fname, fbas, work_dir=None, ignoreDistance=[], with_protein_cost=False, constraint_mode='strict', moma=False):
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
        assert numpy.alltrue((fC == fC[0])), '\nModels must have the same number of fluxes\n!'

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
    #ignoreDistance = []

    Combi = CBTools.ComboGen()
    Cnumber = 2
    Cdata = ''
    unique_combinations = None
    for x in range(len(fbas)):
        Cdata += '%s' % x
    if __DEBUG__: print(Cdata)

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
        if __DEBUG__: print(uq)
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
            if __DEBUG__: print(RiD1[s], RiD2[s])
            if RiD1[s] not in ignoreDistance:
                if moma:
                    ridx2 = RiD2.index(MD1.prefix+RiD1[s])
                else:
                    ridx2 = s
                ##  av = 'z%s' % (zbase+s+1)
                av = 'zvar%s%s' % (combcount,RiD1[s].replace(MD0.prefix,''))

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
    assert len(f.objectives[f.activeObjIdx].getFluxObjectiveReactions()) == 1, "\nOnly single fluxObjectives dealt with at this time"
    if moma:
        f = fbas[0]
        objFcnstr.append('C_%s: %s >= %f' % (f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],\
                                             f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0], f.objectives[f.activeObjIdx].value))
    else:
        for f in fbas:
            objFcnstr.append('C_%s: %s >= %f' % (f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],\
                                                 f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0], f.objectives[f.activeObjIdx].value))
    conL = conL + objFcnstr

    if __DEBUG__:
        print(objL)
        for c in conL:
            print(c)
        print(' ')
        for b in bndL:
            print(b)

    F = open(os.path.join(work_dir, fname+'.lp'), 'w')
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
    return os.path.join(work_dir, fname+'.lp')

def writeMinDistanceLP_absL1(fname, fbas, work_dir=None, ignoreDistance=[], bigM=500, with_protein_cost=False, constraint_mode='strict', moma=False):
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
        assert numpy.alltrue((fC == fC[0])), '\nModels must have the same number of fluxes\n!'


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
    if __DEBUG__: print(Cdata)

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
    ILPMETHOD = 'GK' # gunnar
    zbase = 0
    combcount = 1

    artVar = []
    artVarX = []
    #ignoreDistance = []
    boolVars = []

    for uq in unique_combinations:
        if __DEBUG__: print(uq)
        MD0 = fbas[uq[0]]
        if moma:
            MD1 = fbas[uq[1]].clone()
            MD1.setPrefix('mt_', 'all')
        else:
            MD1 = fbas[uq[1]]

        RiD1 = MD0.getReactionIds()
        RiD2 = MD1.getReactionIds()

        for s in range(len(MD0.reactions)):
            if __DEBUG__: print(RiD1[s], RiD2[s])
            if RiD1[s] not in ignoreDistance:
                if moma:
                    ridx2 = RiD2.index(MD1.prefix+RiD1[s])
                else:
                    ridx2 = s
                Var1 = RiD1[s]
                bVar1 =  'xvar_%s' % Var1
                absVar1 = 'absL_%s' % Var1
                Var2 = RiD2[ridx2]
                bVar2 =  'xvar_%s' % Var2
                absVar2 = 'absL_%s' % Var2

                c0a = '\n'
                c0a += '%s - %s >= 0\n' % (absVar1, Var1)
                c0a += '%s + %s >= 0\n' % (absVar1, Var1)

                if ILPMETHOD == 'GK':
                    # gunnar
                    c0a += '%s + %s - %s %s <= 0\n' % (absVar1, Var1, bigM, bVar1)
                    c0a += '%s - %s + %s %s <= %s\n' % (absVar1, Var1, bigM, bVar1, bigM)
                elif ILPMETHOD == 'SK':
                    # steven
                    c0a += '%s - %s - %s %s <= 0\n' % (absVar1, Var1, bigM, bVar1)
                    c0a += '%s + %s + %s %s <= %s\n' % (absVar1, Var1, bigM, bVar1, bigM)

                c0b = '\n'
                c0b += '%s - %s >= 0\n' % (absVar2, Var2)
                c0b += '%s + %s >= 0\n' % (absVar2, Var2)
                if ILPMETHOD == 'GK':
                    # gunnar
                    c0b += '%s + %s - %s %s <= 0\n' % (absVar2, Var2, bigM, bVar2)
                    c0b += '%s - %s + %s %s <= %s\n' % (absVar2, Var2, bigM, bVar2, bigM)
                elif ILPMETHOD == 'SK':
                    # steven
                    c0b += '%s - %s - %s %s <= 0\n' % (absVar2, Var2, bigM, bVar2)
                    c0b += '%s + %s + %s %s <= %s\n' % (absVar2, Var2, bigM, bVar2, bigM)

                c0 = c0a + c0b

                av = 'zvar%s%s' % (combcount,Var1.replace(MD0.prefix,''))

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
    assert len(f.objectives[f.activeObjIdx].getFluxObjectiveReactions()) == 1, "\nOnly single fluxObjectives dealt with at this time"
    if moma:
        #OFvalue = f.objectives[f.activeObjIdx].value
        f = fbas[0]
        objFcnstr.append('C_%s: %s >= %f' % (f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],\
                                             f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0], f.objectives[f.activeObjIdx].value))
    else:
        for f in fbas:
            objFcnstr.append('C_%s: %s >= %f' % (f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],\
                                                 f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0], f.objectives[f.activeObjIdx].value))
    conL = conL + objFcnstr
    if __DEBUG__:
        print(objL)
        for c in conL:
            print(c)
        print(' ')
        for b in bndL:
            print(b)

    F = open(os.path.join(work_dir, fname+'.lp'), 'w')
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
    return os.path.join(work_dir, fname+'.lp')



def writeMinDistanceLPwithCost(fname, fbas, work_dir=None, ignoreDistance=[], constraint_mode='strict'):
    """
    For backwards compatability only
    """
    print("\n\n**********\nDeprecation warning!\nPlease use writeMinDistanceLP(with_protein_cost=True) instead of writeMinDistanceLPwithCost()\n\n**********\n")
    writeMinDistanceLP(fname, fbas, work_dir=work_dir, ignoreDistance=ignoreDistance, with_protein_cost=True, constraint_mode=constraint_mode)

# CAN GO SOON JUST PUTTING IN FOR SVN SYNCH
"""
def writeMinDistanceLPwithCost(fname, fbas, work_dir=None, ignoreDistance=[], constraint_mode='strict'):
    if work_dir == None:
        work_dir = os.getcwd()

    fC = []
    objFname = ''
    for l in fbas:
        fC.append(len(l.reactions))
        objFname += l.prefix
    fC = numpy.array(fC)
    if __DEBUG__:
        print fC
        print (fC == fC[0])
        print numpy.alltrue((fC == fC[0]))
    assert numpy.alltrue((fC == fC[0])), '\nModels must have the same number of fluxes\n!'

    conL = []
    # model flux_bounds
    initial_cnstr = []

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
    ignoreDistance = []

    Combi = CBTools.ComboGen()
    Cnumber = 2
    Cdata = ''
    unique_combinations = None
    for x in range(len(fbas)):
        Cdata += '%s' % x
    if __DEBUG__: print Cdata

    Combi.uniqueCombinations(Cdata, Cnumber, temp=[])
    Combi.numberifyComb2Int()
    unique_combinations = Combi.combo_int

    if __DEBUG__:
        print 'Data (%s):\n%s\n' % (Cnumber, Cdata)
        print 'UniqueCombinStr:\n%s' % Combi.combo
        print 'UniqueCombinations:\n%s' % unique_combinations

    zbase = 0
    combcount = 1
    for uq in unique_combinations:
        if __DEBUG__: print uq
        RiD1 = fbas[uq[0]].getReactionIds()
        RiD2 = fbas[uq[1]].getReactionIds()
        for s in range(len(fbas[uq[0]].reactions)):
            if __DEBUG__: print RiD1[s], RiD2[s]
            if RiD1[s] not in ignoreDistance:
                ##  av = 'z%s' % (zbase+s+1)
                ##  av = 'zvar%s' % (RiD1[s].replace(fbas[uq[0]].prefix,''))
                av = 'zvar%s%s' % (combcount,RiD1[s].replace(fbas[uq[0]].prefix,''))
                c1 = '%sa: %s - %s - %s <= 0.0' % (av, RiD1[s], RiD2[s], av)
                ##  c1 = '%s - %s - %s <= 0.0' % (RiD1[s], RiD2[s], av)
                c2 = '%sb: %s - %s + %s >= 0.0' % (av, RiD1[s], RiD2[s], av)
                ##  c2 = '%s - %s + %s >= 0.0' % (RiD1[s], RiD2[s], av)
                conL.append(c1)
                conL.append(c2)
                ##  print RiD1[s], RiD2[s], av
                ##  print fbas[uq[0]].reactions[s].getId(), fbas[uq[1]].reactions[s].getId()
                ##  print fbas[uq[0]].reactions[s].annotation['CBM_PEPTIDE_COST'], fbas[uq[1]].reactions[s].annotation['CBM_PEPTIDE_COST']
                av = '%s %s' % (fbas[uq[0]].reactions[s].annotation['CBM_PEPTIDE_COST'], av)
                ##  print RiD1[s], RiD2[s], av
                artVar.append(av)
            else:
                # THIS MAY BE A BUG
                ##  av = 'z%s' % (zbase+s+1)
                ##  av = 'zvar%s' % (RiD1[s].replace(fbas[uq[0]].prefix,''))
                av = 'zvar%s%s' % (combcount,RiD1[s].replace(fbas[uq[0]].prefix,''))
                artVarX.append(av)
                c1 = '\\* %sa: %s - %s - %s <= 0.0 *\\' % (av, RiD1[s], RiD2[s], av)
                ##  c1 = '\\* %s - %s - %s <= 0.0 *\\' % (RiD1[s], RiD2[s], av)
                c2 = '\\* %sb: %s - %s + %s >= 0.0 *\\' % (av, RiD1[s], RiD2[s], av)
                ##  c2 = '\\* %s - %s + %s >= 0.0 *\\' % (RiD1[s], RiD2[s], av)
                conL.append(c1)
                conL.append(c2)
        combcount += 1
        zbase += len(fbas[uq[0]].reactions)

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
    for f in fbas:
        objFcnstr.append('C_%s: %s >= %f' % (f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0],\
                                    f.objectives[f.activeObjIdx].getFluxObjectiveReactions()[0], f.objectives[f.activeObjIdx].value))
    conL = conL + objFcnstr

    if __DEBUG__:
        print objL
        for c in conL:
            print c
        print ' '
        for b in bndL:
            print b

    F = open(os.path.join(work_dir, fname+'.lp'), 'w')
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
    print 'LP written to: %s.lp' % os.path.join(work_dir, fname)
    return os.path.join(work_dir, fname+'.lp')
"""

def writeOptimalSolution(fba, fname, Dir=None, separator=',', only_exchange=False):
    """
    This function writes the optimal solution to file

     - *fba* an instance of an PySCeSCBM model
     - *fname* the output filename
     - *Dir* [default=None] use current directory if not None
     - *separator* [default=','] the column separator
     - *only_exchange* [default=False] only output fluxes labelled as exchange reactions

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
        print('\nCSV file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(fname_r))
        return
    cntr = 0
    F.write('%s%s%s%s%s%s\n' % ('ObjectiveFunction',separator,objName,separator,separator,separator))
    F.write('\"%s\"%s%s%s%s%s%s%s%s%s\"%s\"%s\"%s\"\n' % ('Reaction',separator,'Value',separator,'LowerBound',separator,'UpperBound',separator,'Reduced cost',separator,'Name',separator,'Gene association'))
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
            F.write('\"%s\"%s%s%s%s%s%s%s%s%s\"%s\"%s\"%s\"\n' % (r.getId(),separator,r.value,separator,Lbnd,separator,Ubnd,separator,r.reduced_cost,separator,r.name,separator,gene))
    F.flush()
    F.close()
    print('Reactions exported to {}'.format(fname_r))

def writeModelInfoToFile(fba, fname, Dir=None, separator=',', only_exchange=False, met_type='all'):
    """
    This function writes a CBModel to file

     - *fba* an instance of an PySCeSCBM model
     - *fname* the output filename
     - *Dir* [default=None] use directory if not None
     - *separator* [default=','] the column separator
     - *only_exchange* [default=False] only output fluxes labelled as exchange reactions
     - *type* [default='all'] only output certain type of species: 'all','boundary' or 'variable'

    """
    writeReactionInfoToFile(fba, fname, Dir=Dir, separator=separator, only_exchange=only_exchange)
    writeSpeciesInfoToFile(fba, fname, Dir=Dir, separator=separator, met_type=met_type)

def writeReactionInfoToFile(fba, fname, Dir=None, separator=',', only_exchange=False):
    """
    This function writes a CBModel to file

     - *fba* an instance of an PySCeSCBM model
     - *fname* the output filename
     - *Dir* [default=None] use directory if not None
     - *separator* [default=','] the column separator
     - *only_exchange* [default=False] only output fluxes labelled as exchange reactions

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
        print('\nOutput file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(fname_r))
        return
    cntr = 0
    F.write('%s%s%s%s%s%s\n' % ('ObjectiveFunction',separator,objName,separator,separator,separator))
    F.write('\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n' % ('Reaction',separator,'LowerBound',separator,'UpperBound',separator,'Name',separator,'Equation',separator,'Gene association'))
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
                equation = ' %s '  % __CBCONFIG__['IRREVERSIBLE_SYMBOL']
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
            F.write('\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n' % (r.getId(),separator,Lbnd,separator,Ubnd,separator,r.name,separator,equation,separator,gene))
            GO = False
    F.flush()
    F.close()
    print('Reactions exported to {}'.format(fname_r))

def writeSpeciesInfoToFile(fba, fname, Dir=None, separator=',', met_type='all'):
    """
    This function writes a CBModel to file

     - *fba* an instance of an PySCeSCBM model
     - *fname* the output filename
     - *Dir* [default=None] use directory if not None
     - *separator* [default=','] the column separator
     - *met_type* [default='all'] only output certain type of species: 'all','boundary' or 'variable'

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
    if met_type not in ['all','boundary','variable']:
        met_type = 'all'

    try:
        F2 = open(fname_s, 'w')
        F3 = open(fname_sr, 'w')
    except IOError:
        print('\nOutput file \"{}\" is locked by an external application (probably Excel) please close file and try again (or use a different filename).'.format(fname_s))
        return

    maxReagentOf = 0
    F2.write('\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n' % ('Species',separator,'Name',separator,'Compartment',separator,'Fixed',separator,'ChemicalFormula',separator,'Charge'))
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
            F2.write('\"%s\"%s%s%s%s%s\"%s\"%s\"%s\"%s\"%s\"\n' % (s.getId(),separator,s.getName(),separator,s.compartment,separator,s.is_boundary,separator,s.chemFormula,separator,s.charge))
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
    Prints the FBA optimal solution to the screen.

     - *fba* an FBA model object
     - *include_all* include all variables

    """
    OFflux = fba.objectives[fba.activeObjIdx].fluxObjectives[0].reaction
    OFvalue = fba.objectives[fba.activeObjIdx].value
    OFSense = fba.objectives[fba.activeObjIdx].operation
    print('\n\n**********\nModel: {}\n\n'.format(fba.getId()))
    print('{} objective: {}\nOptimal value: {}\n\n'.format(OFSense,OFflux,OFvalue))
    if include_all:
        for J in fba.reactions:
            print('{}: {}'.format(J.getId(), J.value))
    print('**********\n')

def exportModel(fba, fname=None, fmt='lp', work_dir=None, use_rational='both'):
    """
    Export the FBA model in different formats:

     - *fba* the FBA model
     - *fname* [default=None] the exported filename if None then `fba.getId()` is used
     - *fmt* [default='lp'] the export format can be one of: 'lp' (CPLEX), 'hformat' (Polyhedra), 'all' (both)
     - *use_rational* [default='both'] if *all* or *hformat* is specified should hformat files be written using rational math or not. The default *both* is the legacy behaviour and writes both.

    Note that 'hformat' ignores 'fname' and only uses fba.getId() this is a legacy behaviour

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
            writeModelHFormatFBA2(fba, fname=fname, work_dir=work_dir, use_rational=True)
        elif use_rational == True:
            writeModelHFormatFBA2(fba, fname=fname, work_dir=work_dir, use_rational=True)
        elif use_rational == False:
            writeModelHFormatFBA2(fba, fname=fname, work_dir=work_dir)

def writeProteinCostToCSV(fba, fname):
    """
    Writes the protein costs 'CBM_PEPTIDE_COST' annotation toa csv file.

     - *fba* an instantiated FBA object
     - *fname* the exported file name

    """
    F = open(fname+'.costs.csv','w')
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
        F.write('%s,%s,%s,%s,%s\n' % (rid,minL,maxL,avg_l,pcost))
    F.flush()
    F.close()
    print('Protein costs written to file: {}'.format(fname+'.csv'))

def WriteFVAtoCSV(id, fva, names, Dir=None, fbaObj=None):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV()\"\n')
    time.sleep(1)
    writeFVAtoCSV(fva, names, id, Dir, fbaObj)

#def writeFVAtoCSV(id, fva, names, Dir=None, fbaObj=None):
    #"""
    #INFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV\"
    #"""
    #print('\nINFO: this method will be deprecated please update your scripts to use \"writeFVAtoCSV\"\n')
    #time.sleep(1)
    #writeFVAtoCSV(id, fva, names, Dir, fbaObj)

def writeFVAtoCSV(fvadata, names, fname, Dir=None, fbaObj=None):
    """
    Takes the resuls of a FluxVariabilityAnalysis method and writes it to a nice
    csv file. Note this method replaces the glpk/cplx_WriteFVAtoCSV methods.

     - *fvadata* FluxVariabilityAnalysis() OUTPUT_ARRAY
     - *names* FluxVariabilityAnalysis() OUTPUT_NAMES
     - *fname* filename_base for the CSV output
     - *Dir* [default=None] if set the output directory for the csv files
     - *fbaObj* [default=None] if supplied adds extra model information into the output tables

    """
    if Dir != None:
        Dir = os.path.join(Dir, fname+'.fva.csv')
    else:
        Dir = fname+'.fva.csv'
    F = open(Dir, 'w')
    if fbaObj == None:
        F.write('name,optval,min,max,diff,red cost,minstat,maxstat\n')
    else:
        F.write('name,optval,min,max,diff,red cost,minstat,maxstat,"equation","subsystem","gene association","confidence level"\n')
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
                F.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (name, optval, min, max, diff, rc, minstat, maxstat))
            else:
                xInf = []
                Ro = fbaObj.getReaction(name)
                for k in ['Equation','SUBSYSTEM','GENE ASSOCIATION','Confidence Level']:
                    if k in Ro.annotation:
                        xInf.append(Ro.annotation[k])
                    else:
                        xInf.append('')
                F.write('%s,%s,%s,%s,%s,%s,%s,%s,"%s","%s","%s",%s\n' % (name, optval, min, max, diff, rc, minstat, maxstat, xInf[0], xInf[1], xInf[2], xInf[3]))

    F.flush()
    F.close()
    print('FVA results written to: {}'.format(Dir))

def WriteFVAdata(fva, names, fname, work_dir=None, roundec=None, scale_min=False, appendfile=False, info=None):
    """
    INFO: this method will be deprecated please update your scripts to use \"writeFVAdata()\"
    """
    print('\nINFO: this method will be deprecated please update your scripts to use \"writeFVAdata()\"\n')
    time.sleep(1)
    writeFVAdata(fva, names, fname, work_dir, roundec, scale_min, appendfile, info)

def writeFVAdata(fvadata, names, fname, work_dir=None, roundec=None, scale_min=False, appendfile=False, info=None):
    """
    Takes the resuls of a FluxVariabilityAnalysis method and writes it to a nice
    csv file. Note this method replaces the glpk/cplx_WriteFVAtoCSV methods. Data is output as a csv file
    with columns: FluxName, FVA_MIN, FVA_MAX, OPT_VAL, SPAN

     - *fvadata* FluxVariabilityAnalysis() FVA OUTPUT_ARRAY
     - *names* FluxVariabilityAnalysis() FVA OUTPUT_NAMES
     - *fname* filename_base for the CSV output
     - *work_dir* [default=None] if set the output directory for the csv files
     - *roundec* [default=None] an integer indicating at which decimal to round off output. Default is no rounding.
     - *scale_min* [default=False] normalise each flux such that that FVA_MIN = 0.0
     - *appendfile* [default=False] instead of opening a new file try and append the data
     - *info* [default=None] a string added to the results as an extra column, useful with `appendfile`

    """
    if work_dir != None:
        work_dir = os.path.join(work_dir, fname+'.fvadata.csv')
    else:
        work_dir = fname+'.fvadata.csv'
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
                    max = max-min
                    optval = optval-min
                    min = 0.0
                elif min <= 0.0:
                    max = max + abs(min)
                    optval = optval + abs(min)
                    min = 0.0
            if info == None:
                if roundec == None:
                    F.write('%s,%s,%s,%s,%s\n' % (name, min, max, optval, abs(fvadata[Jidx][3] - fvadata[Jidx][2])))
                else:
                    F.write('%s,%s,%s,%s,%s\n' % (name, min, max, optval, round(abs(fvadata[Jidx][3] - fvadata[Jidx][2]), roundec)))

            else:
                if roundec == None:
                    F.write('%s,%s,%s,%s,%s,%s\n' % (name, min, max, optval, abs(fvadata[Jidx][3] - fvadata[Jidx][2]), info))
                else:
                    F.write('%s,%s,%s,%s,%s,%s\n' % (name, min, max, optval, round(abs(fvadata[Jidx][3] - fvadata[Jidx][2]), roundec), info))
    F.flush()
    F.close()
    print('FVAdata results written to: {}'.format(work_dir))

def writeSolutions(fname, sols=[], sep=',', extra_output=None, fba=None):
    """
    Write 2 or more solutions where a solution is a dictionary of flux:value pairs:

     - *fname* the export filename
     - *sols* a list of dictionaries containing flux:value pairs (e.g. output by cmod.getReactionValues())
     - *sep* [default=','] the column separator
     - *extra_output* [default=None] add detailed information to output e.g. reaction names by giving a CBModel object as an argument to *extra_output*.
     - *fba* an fba model that canbe used for extra_output

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

    F = open(fname+'.csv', 'w')
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
            row = row[:-1]+'\n'
        F.write(row)
    F.flush()
    F.close()
    print('\nSolutions written to: \"%s\"\n' % fname)


def generateBGID(num, prefix):
    """
    Create a BGID generator, which is <prefix><num> where perfix is two letters num is padded to 6 figures
     - *num* the starting number
     - *prefix* the two letter prefix

    """
    assert prefix[0].isalpha() and prefix[1].isalpha(), '\nInvalid prefix'
    while True:
        num += 1
        yield '{}{:0>6}'.format(prefix, num)

def convertFloatToExcel(num, roundoff):
    """
    Converts a float to Excel compatible "number"

     - *num* a number
     - *roundoff* the number of roundoff digits for round()

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
    Converts an Excel "number" to a float

     - *num* a number

    """
    ## keep this in sync with convertFloatToExcel
    if num == '':
        rval = None
    elif num is 'nan':
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
    Exports the model as an Excel 97 spreadsheet

     - *fba* a CBMPy model instance
     - *filename* the filename of the workbook
     - *roundoff* [default=6] the number of digits to round off to

    """

    if not _HAVE_XLWT_:
        print('\nERROR: Cannot create Excel file, XLWT package not available (http://pypi.python.org/pypi/xlwt)')
        return
    try:
        F = open('{}.xls'.format(filename), 'wb')
        F.close()
    except:
        print('\nERROR: cannot open file "{}"! Please close workbook before writing!\n'.format('{}.xls'.format(filename)))
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
                wsInf.write(ridx, cidx+1, t_[1])
            except:
                wsInf.write(ridx, cidx+1, 'undefined')
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
            wsInf.write(ridx, t_+1, fba.DATE_CREATED[t_])
    ridx += 1
    wsInf.write(ridx, 0, 'lastmodified', styleBold)
    if fba.DATE_MODIFIED != None:
        for t_ in range(len(fba.DATE_MODIFIED)):
            wsInf.write(ridx, t_+1, fba.DATE_MODIFIED[t_])
    ridx += 1
    cTime = list(time.gmtime())
    wsInf.write(ridx, 0, 'exported', styleBold)
    for t_ in range(len(cTime)):
        if t_<= 5:
            wsInf.write(ridx, t_+1, cTime[t_])
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
    wsInf.write(ridx, 0, xlwt.Formula('HYPERLINK("{}";"{}")'.format('http://cbmpy.sourceforge.net', 'Exported by CBMPy {}'.format(__version__))), styleHyper)
    ridx += 1

    #Scol = ['metabolite', 'value', 'shadow price', '', 'reaction', 'flux', 'reduced cost', 'FVA min', 'FVA max']
    if fba.SCALED_REDUCED_COSTS:
        rcstr = 'scaled reduced cost'
    else:
        rcstr = 'reduced cost'

    Scol = ['reaction', 'flux', 'lower', 'upper', rcstr, 'FVA min', 'FVA max', 'FVA span', 'exchange', 'info', 'stoichiometry']
    for s_ in range(len(Scol)):
        wsSol.write(0, s_, Scol[s_], styleBoldC)

    # metabolites
    Mlist = []
    Mcols = ['id','name','charge','chemformula','compartment','fixed']
    #Mcols = ['id','name','charge','chemformula','compartment','fixed','bgid']
    MUcols = []

    MiriQual = ["is","isEncodedBy","encodes","hasPart","hasProperty","hasTaxon","hasVersion","isDescribedBy",\
                "isHomologTo","isPartOf","isPropertyOf","isVersionOf","occursIn"]
    wsMiriam.write(0, 0, 'id', styleBoldC)
    for q_ in range(len(MiriQual)):
        wsMiriam.write(0, q_+1, MiriQual[q_], styleBold)

    # get data
    mcntr = 1
    for s_ in range(len(fba.species)):
        Mdi = {}
        Mdi['id'] = fba.species[s_].getId()
        Mdi['name'] = fba.species[s_].getName()
        Mdi['charge'] = fba.species[s_].getCharge()
        Mdi['chemformula'] = fba.species[s_].getChemFormula()
        Mdi['compartment'] = fba.species[s_].compartment
        Mdi['fixed'] = fba.species[s_].is_boundary
        #Mdi['bgid'] = bgGen.next()
        udata = {}
        for k_ in fba.species[s_].getAnnotations():
            if k_ not in MUcols:
                MUcols.append(k_)
            udata[k_] = fba.species[s_].annotation[k_]
        Mdi['data'] = udata
        Mdi['reagentof'] = fba.species[s_].isReagentOf()
        Mlist.append(Mdi)

        ## TODO: removed for now will appear in its own sheet at some point
        #wsSol.write(s_+1, 0,  Mdi['id'])
        #Sval = fba.species[s_].getValue()
        #if Sval == None or numpy.isnan(Sval) or numpy.isinf(Sval):
            #Sval = ''
        #wsSol.write(s_+1, 1,  Sval)
        #wsSol.write(s_+1, 2,  fba.species[s_].shadow_price)
        annot = fba.species[s_].getMIRIAMannotations()
        if annot != None:
            for k_ in annot:
                for m_ in annot[k_]:
                    wsMiriam.write(mcntr, 0, Mdi['id'])
                    #if 'identifiers.org' in m_:
                        #wsMiriam.write(mcntr, MiriQual.index(k_)+1, xlwt.Formula('HYPERLINK("{}";"{}")'.format(m_.replace('identifiers.org', 'info.identifiers.org'),\
                                                                                                #m_.rsplit('/',1)[1])), styleHyper)
                    #else:
                        #wsMiriam.write(mcntr, MiriQual.index(k_)+1, m_)
                    wsMiriam.write(mcntr, MiriQual.index(k_)+1, m_)
                    mcntr += 1

    # reactions
    Rlist = []
    #Rcols = ['id','name','reversible','lowerbound','upperbound','compartment','bgid']
    Rcols = ['id','name','reversible','lowerbound','upperbound','compartment']
    RUcols = []
    try:
        AOIDS = fba.getActiveObjective().getFluxObjectiveReactions()
    except:
        AOIDS = []
    rids = fba.getReactionIds()
    rids.sort()
    for r_ in range(len(fba.reactions)):
        Rdi = {}
        REAC = fba.getReaction(rids[r_])
        Rdi['id'] = REAC.getId()
        Rdi['name'] = REAC.getName()
        Rdi['reversible'] = REAC.reversible
        Rdi['compartment'] = REAC.compartment
        #Rdi['bgid'] = bgGen.next()
        Rdi['lowerbound'] = fba.getReactionLowerBound(Rdi['id'])
        Rdi['upperbound'] = fba.getReactionUpperBound(Rdi['id'])
        udata = {}
        for k_ in REAC.getAnnotations():
            if REAC.annotation[k_] is not None:
                if k_ not in RUcols:
                    RUcols.append(k_)
                udata[k_] = REAC.annotation[k_]
        Rdi['data'] = udata
        subs = []
        prods = []
        for rr_ in REAC.reagents:
            if rr_.getCoefficient() > 0.0:
                prods.append((abs(rr_.getCoefficient()), rr_.getSpecies()))
            else:
                subs.append((abs(rr_.getCoefficient()), rr_.getSpecies()))
        Rdi['substrates'] = subs
        Rdi['products'] = prods
        Rlist.append(Rdi)
        cstart = 0
        bnds = list(fba.getReactionBounds(Rdi['id']))
        # deal with equality
        if bnds[3] != None:
            bnds[1] = bnds[3]
            bnds[2] = bnds[3]
        if Rdi['id'] not in AOIDS:
            wsSol.write(r_+1, cstart,  Rdi['id'])
        else:
            wsSol.write(r_+1, cstart,  Rdi['id'], styleBold)
        rval = REAC.getValue()
        rval = convertFloatToExcel(rval, roundoff)
        wsSol.write(r_+1, cstart+1, rval)
        del rval

        wsSol.write(r_+1, cstart+2, convertFloatToExcel(bnds[1], roundoff))
        wsSol.write(r_+1, cstart+3, convertFloatToExcel(bnds[2], roundoff))
        wsSol.write(r_+1, cstart+4, convertFloatToExcel(REAC.reduced_cost, roundoff))
        if REAC.fva_min != None:
            wsSol.write(r_+1, cstart+5, convertFloatToExcel(REAC.fva_min, roundoff))
        if REAC.fva_max != None:
            wsSol.write(r_+1, cstart+6, convertFloatToExcel(REAC.fva_max, roundoff))
        if REAC.fva_min != None and REAC.fva_max != None:
            wsSol.write(r_+1, cstart+7, convertFloatToExcel(REAC.fva_max-REAC.fva_min, roundoff))
        if REAC.is_exchange:
            wsSol.write(r_+1, cstart+8, 'yes')
        #else:
            #wsSol.write(r_+1, cstart+8, 'no')
        wsSol.write(r_+1, cstart+9, xlwt.Formula('HYPERLINK("#reactions!A{}";"{}")'.format(r_+1, 'info')), styleHyper)
        wsSol.write(r_+1, cstart+10, xlwt.Formula('HYPERLINK("#network_react!A{}";"{}")'.format(((r_+1)*3)-2, 'stoich')), styleHyper)

        annot = REAC.getMIRIAMannotations()
        if annot != None:
            for k_ in annot:
                for m_ in annot[k_]:
                    wsMiriam.write(mcntr, 0, Rdi['id'])
                    #if 'identifiers.org' in m_:
                        #wsMiriam.write(mcntr, MiriQual.index(k_)+1, xlwt.Formula('HYPERLINK("{}";"{}")'.format(m_.replace('identifiers.org', 'info.identifiers.org'),\
                                                                                                #m_.rsplit('/',1)[1])), styleHyper)
                    #else:
                        #wsMiriam.write(mcntr, MiriQual.index(k_)+1, m_)
                    wsMiriam.write(mcntr, MiriQual.index(k_)+1, m_)
                    mcntr += 1


    compHead = ['id', 'name', 'size', 'dimensions', '# species', '# reactions']
    compHeadU = []
    for c_ in range(len(compHead)):
        wsComp.write(0, c_, compHead[c_], styleBoldC)
    ridx = 1
    for c_ in fba.compartments:
        wsComp.write(ridx, 0, c_.getId())
        wsComp.write(ridx, 1, c_.getName())
        wsComp.write(ridx, 2, convertFloatToExcel(c_.getSize(), roundoff))

        wsComp.write(ridx, 3, c_.getDimensions())
        wsComp.write(ridx, 4, len(c_.containsSpecies()))
        cRe = c_.containsReactions()
        if len(cRe) != 0:
            wsComp.write(ridx, 5, len(cRe))

        annot = c_.getAnnotations()
        for k_ in annot:
            if k_ not in compHeadU:
                compHeadU.append(k_)
        for k_ in annot:
            wsComp.write(ridx, 6+compHeadU.index(k_), annot[k_])

        annot = c_.getMIRIAMannotations()
        if annot != None:
            for k_ in annot:
                for m_ in annot[k_]:
                    wsMiriam.write(mcntr, 0, c_.getId())
                    wsMiriam.write(mcntr, MiriQual.index(k_)+1, m_)
                    mcntr += 1
        ridx += 1
    for ud_ in range(0, len(compHeadU)):
        wsComp.write(0, 6+ud_, compHeadU[ud_], styleItalicC)
    ridx = 1

    del ridx, c_

    # write data
    # write metabolites
    for m_ in range(len(Mcols)):
        wsMet.write(0, m_, Mcols[m_], styleBold)
    for m_ in range(len(MUcols)):
        wsMet.write(0, len(Mcols)+m_, MUcols[m_], styleItalic)
    stridx = 0
    for s_ in range(len(Mlist)):
        s = s_+1
        wsMet.write(s, Mcols.index('id'), Mlist[s_]['id'])
        wsMet.write(s, Mcols.index('name'), Mlist[s_]['name'])
        wsMet.write(s, Mcols.index('charge'), Mlist[s_]['charge'])
        wsMet.write(s, Mcols.index('chemformula'), Mlist[s_]['chemformula'])
        wsMet.write(s, Mcols.index('compartment'), Mlist[s_]['compartment'])
        wsMet.write(s, Mcols.index('fixed'), Mlist[s_]['fixed'])
        #wsMet.write(s, Mcols.index('bgid'), Mlist[s_]['bgid'])
        for ud_ in Mlist[s_]['data']:
            try:
                if Mlist[s_]['data'][ud_] is not None:
                    if len(Mlist[s_]['data'][ud_]) < 30000:
                        wsMet.write(s, len(Mcols)+MUcols.index(ud_), str(Mlist[s_]['data'][ud_]))
                    else:
                        wsMet.write(s, len(Mcols)+MUcols.index(ud_),'Data too long (more than 30000 characters)', styleBold)
            except TypeError:
                print('Annotation write error (TypeError) {} : {}'.format(ud_, Mlist[s_]['data'][ud_]))

        wsStR.write(stridx, 0, Mlist[s_]['id'], styleBold)
        cidx = 0
        stridx += 1
        for ro_ in range(len(Mlist[s_]['reagentof'])):
            wsStR.write(stridx, cidx, Mlist[s_]['reagentof'][ro_])
            cidx += 1
            if cidx >= 200:
                cidx = 0
                stridx += 1
        stridx += 1

    # write reactions
    for r_ in range(len(Rcols)):
        wsRe.write(0, r_, Rcols[r_], styleBold)
    for r_ in range(len(RUcols)):
        try:
            wsRe.write(0, len(Rcols)+r_, RUcols[r_], styleItalic)
        except Exception as why:
            print(why)
            print('\nERROR: Exceeded Excel columns')
    ridx = 0
    for r_ in range(len(Rlist)):
        r = r_+1
        wsRe.write(r, Rcols.index('id'), Rlist[r_]['id'])
        wsRe.write(r, Rcols.index('name'), Rlist[r_]['name'])
        wsRe.write(r, Rcols.index('reversible'), Rlist[r_]['reversible'])
        wsRe.write(r, Rcols.index('lowerbound'), convertFloatToExcel(Rlist[r_]['lowerbound'], roundoff))
        wsRe.write(r, Rcols.index('upperbound'), convertFloatToExcel(Rlist[r_]['upperbound'], roundoff))
        wsRe.write(r, Rcols.index('compartment'), Rlist[r_]['compartment'])
        #wsRe.write(r, Rcols.index('bgid'), Rlist[r_]['bgid'])
        for ud_ in Rlist[r_]['data']:
            if Rlist[r_]['data'][ud_] == None:
                wsRe.write(r, len(Rcols)+RUcols.index(ud_), '')
            elif type(Rlist[r_]['data'][ud_]) == str and len(Rlist[r_]['data'][ud_]) >= 30000:
                wsRe.write(r, len(Rcols)+RUcols.index(ud_), 'Data too long (more than 30000 characters)', styleBold)
            else:
                try:
                    wsRe.write(r, len(Rcols)+RUcols.index(ud_), Rlist[r_]['data'][ud_])
                except:
                    print('Exceeded Excel columns')

        wsSt.write(ridx, 0, Rlist[r_]['id'], styleBold)
        ridx += 1
        wsSt.write(ridx, 0, 'substrates', styleBold)
        cidx = 1
        for sr_ in range(len(Rlist[r_]['substrates'])):
            wsSt.write(ridx, cidx+sr_, Rlist[r_]['substrates'][sr_][0])
            wsSt.write(ridx, cidx+sr_+1, Rlist[r_]['substrates'][sr_][1])
            cidx += 1
        ridx += 1
        wsSt.write(ridx, 0, 'products', styleBold)
        cidx = 1
        for pr_ in range(len(Rlist[r_]['products'])):
            wsSt.write(ridx, cidx+pr_, Rlist[r_]['products'][pr_][0])
            wsSt.write(ridx, cidx+pr_+1, Rlist[r_]['products'][pr_][1])
            cidx += 1
        ridx += 1

    gids = fba.getGroupIds()
    for g_ in range(len(gids)):
        mbrs = fba.getGroup(gids[g_]).getMemberIDs()
        wsGrp.write(0, g_, '{} ({})'.format(gids[g_], fba.getGroup(gids[g_]).getName()), styleBold)
        for m_ in range(len(mbrs)):
            wsGrp.write(m_+1, g_, str(mbrs[m_]))

    wb.save('{}.xls'.format(filename))

def writeModelToCOMBINEarchive(mod, fname=None, directory=None, sbmlname=None, withExcel=True, vc_given='CBMPy', vc_family='Software', vc_email='None', vc_org='cbmpy.sourceforge.net', add_cbmpy_annot=True, add_cobra_annot=True):
    """
    Write a model in SBML and Excel format to a COMBINE archive using the following information:

    - *mod* a model object
    - *fname* the output base filename, archive will be <fname>.zip
    - *directory* [default=None] created the combine archive 'directory'
    - *sbmlname* [default='None'] If *sbmlname* is defined then SBML file is <sbmlname>.xml otherwise sbml will be <fname>.xml.
    - *withExcel* [default=True] include a human readable Excel spreadsheet version of the model
    - *vc_given* [default='CBMPy'] first name
    - *vc_family* [default='Software'] family name
    - *vc_email* [default='None'] email
    - *vc_org* [default='None'] organisation
    - *add_cbmpy_annot* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
    - *add_cobra_annot* [default=True] add COBRA <notes> annotation

    """
    scTime = time.strftime('%Y-%m-%dT%H:%M:%S') + '%i:00' % (time.timezone/60/60)
    #self.writeSedXML(sedx=True)
    #sedxname = '%s.sed.omex' % (self.id)
    #sf = os.path.join(self.sedpath, sedxname)
    #self.__sedarchive__ = sf
    if directory != None:
        zfpath = os.path.join(directory, fname+'.sbex.zip')
    else:
        zfpath = fname+'.sbex.zip'
    zf = zipfile.ZipFile(zfpath, mode='w', compression=zipfile.ZIP_DEFLATED)
    if sbmlname != None:
        if sbmlname.endswith('.xml'):
            sbmlf = sbmlname
            xlf = sbmlname[:-4]
        else:
            sbmlf = sbmlname+'.xml'
            xlf = sbmlname
    else:
        sbmlf = fname+'.xml'
        xlf = fname
    ptmp = os.path.join(os.getcwd(), 'sedxtmp')
    if not os.path.exists(ptmp):
        os.makedirs(ptmp)
    assert os.path.exists(ptmp), "Could not create temporary archive directory: {}".format(ptmp)

    MFstr = ''
    MDstr = ''
    MFstr += '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">\n'
    MFstr += ' <content location="." format="http://identifiers.org/combine.specifications/omex"/>\n'
    #MFstr += ' <content location="./%s" format="http://identifiers.org/combine.specifications/sedml"/>\n' % os.path.split(self.__sedxml__)[-1]
    MFstr += ' <content location="./metadata.rdf" format="http://identifiers.org/combine.specifications/omex-metadata"/>\n'

    # SBML
    writeSBML3FBCV2(mod, sbmlf, ptmp, add_cbmpy_annot=add_cbmpy_annot, add_cobra_annot=add_cobra_annot)
    zf.write(os.path.join(ptmp, sbmlf), arcname=sbmlf)
    MFstr += ' <content location="./{}" format="http://identifiers.org/combine.specifications/sbml.level-3.version-1"/>\n'.format(sbmlf)

    # Excel
    if withExcel and _HAVE_XLWT_:
        writeModelToExcel97(mod, os.path.join(ptmp, xlf))
        xlf += '.xls'
        zf.write(os.path.join(ptmp, xlf), arcname=xlf)
        MFstr += ' <content location="./{}" format="http://mediatypes.appspot.com/application/vnd.ms-excel"/>'.format(xlf)
    MF = open(os.path.join(ptmp, 'manifest.xml'), 'w')
    MF.write('<?xml version="1.0" encoding="utf-8"?>\n{}\n</omexManifest>\n'.format(MFstr))
    MF.close()


    MD = open(os.path.join(ptmp, 'metadata.rdf'), 'w')
    MD.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    MD.write('<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n')
    MD.write('    xmlns:dcterms="http://purl.org/dc/terms/"\n')
    MD.write('    xmlns:vCard="http://www.w3.org/2006/vcard/ns#"\n')
    MD.write('    xmlns:bqmodel="http://biomodels.net/models-qualifiers">\n')
    MD.write(' <rdf:Description rdf:about=".">\n')
    #MDstr += '   <dcterms:description>\n     %s\n    </dcterms:description>\n' % self.omex_description
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
    print('COMBINE archive created: {}'.format(fname+'.zip'))
