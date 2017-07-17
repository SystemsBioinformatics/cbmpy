"""
CBMPy: CBRead module
====================
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
Last edit: $Author: bgoli $ ($Id: CBRead.py 604 2017-07-17 12:30:08Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

#__doc__ = """
#CBMPy: CBRead module
#"""


import os, time, numpy
# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio
from . import CBXML, CBModel

from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

##  try:
    ##  import psyco
    ##  psyco.full()
##  except:
    ##  pass

_HAVE_SYMPY_ = None
try:
    import sympy
    if int(sympy.__version__.split('.')[1]) >= 7 and int(sympy.__version__.split('.')[2]) >= 4:
        HAVE_SYMPY = True
    else:
        del sympy
except ImportError:
    print('Rational IO not available')
    _HAVE_SYMPY_ = False
_HAVE_SYMPY_ = None
try:
    import h5py
    _HAVE_HD5_ = True
except ImportError:
    _HAVE_HD5_ = False

_HAVE_XLRD_ = False
try:
    import xlrd
    _HAVE_XLRD_ = True
except ImportError:
    print('\nINFO: No xlrd module available, Excel spreadsheet reading disabled')



__example_models__ = {'cbmpy_test_core' : 'core_memesa_model.l3.xml',
                   'cbmpy_test_ecoli' : 'Ecoli_iJR904.glc.l3.xml',
                   }
__example_model_path__ = os.path.join(__CBCONFIG__['CBMPY_DIR'], 'models')

if not os.path.exists(os.path.join(__example_model_path__, 'core_memesa_model.l3.xml')) or\
   not os.path.exists(os.path.join(__example_model_path__, 'Ecoli_iJR904.glc.l3.xml')):
    import zipfile
    print('Installing default models ...')
    zfile = zipfile.ZipFile(os.path.join(__example_model_path__, 'default_models.zip.py'), allowZip64=True)
    zfile.extractall(path=__example_model_path__)
    zfile.close()
    del zipfile, zfile


def readSBML3FBC(fname, work_dir=None, return_sbml_model=False, xoptions={'validate' : False}, scan_notes_gpr=True):
    """
    Read in an SBML Level 3 file with FBC annotation where and return a CBM model object

     - *fname* is the filename
     - *work_dir* is the working directory
     - *return_sbml_model* deprecated and ignored please update code
     - *xoptions* special load options, enable with option=True except for nmatrix_type which has a type.

       - *nogenes* do not load/process genes
       - *noannot* do not load/process any annotations
       - *validate* validate model and display errors and warnings before loading
       - *readcobra* read the cobra annotation
       - *read_model_string* [default=False] read the model from a string (instead of a filename) containing an SBML document
       - *nmatrix_type* [default='normal'] define the type of stoichiometrich matrix to be built

         - 'numpy' dense numpy array (best performance)
         - 'scipy_csr' scipy sparse matrix (lower performance, low memory)
         - 'sympy' a sympy rational matrix (low performance, high memory, cast to dense to analyse)
         - None do not build matrix

    - *scan_notes_gpr* [default=True] if the model is loaded and no genes are detected scan the <notes> field for GPR associationa

    """
    if fname in __example_models__:
        fname = __example_models__[fname]
        print(fname)
        fname = os.path.join(__example_model_path__, fname)

    xmod = CBXML.sbml_readSBML3FBC(fname, work_dir, return_sbml_model, xoptions)
    if scan_notes_gpr and len(xmod.getGeneIds()) == 0:
        print('INFO: no standard gene encoding detected, attempting to load from annotations.')
        xmod.createGeneAssociationsFromAnnotations()
    return xmod

def readCOBRASBML(fname, work_dir=None, return_sbml_model=False, delete_intermediate=False, fake_boundary_species_search=False,\
                  output_dir=None, skip_genes=False, scan_notes_gpr=True):
    """
    Read in a COBRA format SBML Level 2 file with FBA annotation where and return either a CBM model object
    or a (cbm_mod, sbml_mod) pair if return_sbml_model=True

     - *fname* is the filename
     - *work_dir* is the working directory
     - *delete_intermediate* [default=False] delete the intermediate SBML Level 3 FBC file
     - *fake_boundary_species_search* [default=False] after looking for the boundary_condition of a species search for overloaded id's <id>_b
     - *output_dir* [default=None] the directory to output the intermediate SBML L3 files (if generated) default to input directory
     - *skip_genes* [default=False] do not load GPR data
     - *scan_notes_gpr* [default=True] if the model is loaded and no genes are detected the scan the <notes> field for GPR associationa

    """
    xmod = CBXML.sbml_readCOBRASBML(fname, work_dir=work_dir, return_sbml_model=False, delete_intermediate=delete_intermediate, fake_boundary_species_search=fake_boundary_species_search, output_dir=output_dir, skip_genes=skip_genes)
    if scan_notes_gpr and len(xmod.getGeneIds()) == 0:
        xmod.createGeneAssociationsFromAnnotations()
    return xmod


def readSBML2FBA(fname, work_dir=None, return_sbml_model=False, fake_boundary_species_search=False, scan_notes_gpr=True):
    """
    Read in an SBML Level 2 file with FBA annotation where:

     - *fname* is the filename
     - *work_dir* is the working directory if None then only fname is used
     - *return_sbml_model* [default=False] return a a (cbm_mod, sbml_mod) pair
     - *fake_boundary_species_search* [default=False] after looking for the boundary_condition of a species search for overloaded id's <id>_b
     - *scan_notes_gpr* [default=True] if the model is loaded and no genes are detected the scan the <notes> field for GPR associationa

    """

    xmod = CBXML.sbml_readSBML2FBA(fname, work_dir, return_sbml_model=False, fake_boundary_species_search=fake_boundary_species_search)
    if scan_notes_gpr and len(xmod.getGeneIds()) == 0:
        xmod.createGeneAssociationsFromAnnotations()
    return xmod


def readLPtoList(fname, work_dir):
    NEW = False
    TYPE = None
    Object = []
    Constr =[]
    Bounds = []
    F = file(os.path.join(work_dir, fname), 'r')
    for l in F:
        if l == '' or l[:2] == '\\\\' or l == '\n' or l.strip() == 'END':
            print('skipping')
            if __DEBUG__: print(l)
        else:
            L = l.strip()
            if L == 'Maximize':
                TYPE = 'ObjFunc'
                NEW = True
            elif L == 'Subject To':
                TYPE = 'Constr'
                NEW = True
            elif L == 'Bounds':
                TYPE = 'Bounds'
                NEW = True
            else:
                NEW = False

            if TYPE == 'ObjFunc' and not NEW:
                Object.append(L)
            elif TYPE == 'Constr' and not NEW:
                Constr.append(L)
            elif TYPE == 'Bounds' and not NEW:
                Bounds.append(L)

    F.close()
    if __DEBUG__:
        print('ObjectiveLines')
        print(Object)
        print('ConstraintLines')
        print(Constr)
        print('BoundsLines')
        print(Bounds)
    return Object, Constr, Bounds

def readSK_FVA(filename):
    """
    Read Stevens FVA results (opt.fva) file and return a list of dictionaries

    """
    assert os.path.exists(filename), '\nGive me a break!\n'
    ##  name = []
    vari = []
    F = file(filename, 'r')
    for l in F:
        L = l.split(':')
        Jn = L[0].strip()
        Vmin = L[1].strip()
        V = L[2].strip()
        V = V.split('--')
        Vmax = V[0].strip()
        Vstat = V[1].strip()
        ##  name.append((Jn, Vstat))
        ##  vari.append((Vmin, Vmax))
        vari.append({'name' : Jn,
                     'min' : Vmin,
                     'max' : Vmax,
                     'status' : Vstat
                     })
        if __DEBUG__: print(Jn, Vmin, Vmax, Vstat)
    return vari

##  def readSK_vertexOld(fname, bigfile=False):
    ##  """
    ##  Reads in Stevens vertex analysis file and returns:

        ##  - a list of vertex vectors
        ##  - a list of ray vectors
        ##  - the basis of the lineality space as a list of vectors

    ##  all vectors in terms of the column space of N

    ##  """

    ##  assert _HAVE_SYMPY_, 'Install Sympy for rational IO support'

    ##  assert os.path.exists(fname), 'Uhm exqueeze me ...'
    ##  SK_vert_file = file(fname, 'r')
    ##  VertOut = []
    ##  LinOut = []
    ##  RayOut = []
    ##  if bigfile:
        ##  VertTmp = gzip.open('_vtx_.tmp.gz','wb', compresslevel=3)
        ##  LinTmp = gzip.open('_lin_.tmp.gz','wb', compresslevel=3)
        ##  RayTmp = gzip.open('_ray_.tmp.gz','wb', compresslevel=3)
    ##  GOvert = False
    ##  GOray = False
    ##  GOlin = False
    ##  lcntr = 0
    ##  lcntrtmp = 0
    ##  for l in SK_vert_file:
        ##  lcntr += 1
        ##  if lcntr == 1000:
            ##  print 'Processing vertex: %s' % (lcntr + lcntrtmp)
            ##  lcntrtmp += lcntr
            ##  lcntr = 0
        ##  if '* Lineality basis ' in l:
            ##  GOvert = False
            ##  GOray = False
            ##  GOlin = True
        ##  if '* Rays ' in l:
            ##  GOvert = False
            ##  GOray = True
            ##  GOlin = False
        ##  if '* Vertices ' in l:
            ##  GOvert = True
            ##  GOray = False
            ##  GOlin = False
        ##  if l[:2] != '* ':
            ##  L = l.split()
            ##  rowL = []
            ##  for c in L:
                ##  c = c.strip()
                ##  if c == '0':
                    ##  rnum = '0'
                    ##  rowL.append(0.0)
                ##  else:
                    ##  rnum = sympy.Rational('%s' % c)
                    ##  rowL.append(rnum.evalf())
                ##  del rnum
            ##  del L
            ##  if GOlin:
                ##  if bigfile:
                    ##  rowEnd = len(rowL)
                    ##  cntr = 0
                    ##  for e in rowL:
                        ##  cntr += 1
                        ##  if e == 0.0:
                            ##  LinTmp.write('0.0')
                        ##  else:
                            ##  LinTmp.write('%.14f' % e)
                        ##  if cntr == rowEnd:
                            ##  LinTmp.write('\n')
                        ##  else:
                            ##  LinTmp.write(',')
                ##  else:
                    ##  LinOut.append(rowL)
            ##  elif GOray:
                ##  if bigfile:
                    ##  rowEnd = len(rowL)
                    ##  cntr = 0
                    ##  for e in rowL:
                        ##  cntr += 1
                        ##  if e == 0.0:
                            ##  RayTmp.write('0.0')
                        ##  else:
                            ##  RayTmp.write('%.14f' % e)
                        ##  if cntr == rowEnd:
                            ##  RayTmp.write('\n')
                        ##  else:
                            ##  RayTmp.write(',')
                ##  else:
                    ##  RayOut.append(rowL)
            ##  elif GOvert:
                ##  if bigfile:
                    ##  rowEnd = len(rowL)
                    ##  cntr = 0
                    ##  for e in rowL:
                        ##  cntr += 1
                        ##  if e == 0.0:
                            ##  VertTmp.write('0.0')
                        ##  else:
                            ##  VertTmp.write('%.14f' % e)
                        ##  if cntr == rowEnd:
                            ##  VertTmp.write('\n')
                        ##  else:
                            ##  VertTmp.write(',')
                ##  else:
                    ##  VertOut.append(rowL)
            ##  del rowL
    ##  print '\nProcessed %s vertices.\n' % (lcntr + lcntrtmp)
    ##  SK_vert_file.close()
    ##  if bigfile:
        ##  VertTmp.close()
        ##  RayTmp.close()
        ##  LinTmp.close()
        ##  VertTmp = gzip.open('_vtx_.tmp.gz','rb')
        ##  LinTmp = gzip.open('_lin_.tmp.gz','rb')
        ##  RayTmp = gzip.open('_ray_.tmp.gz','rb')
        ##  return VertTmp, RayTmp, LinTmp
    ##  else:
        ##  print 'Lineality basis: %s' % len(LinOut)
        ##  print 'Number of rays: %s' % len(RayOut)
        ##  print 'Number of vertices: %s' % len(VertOut)
        ##  return VertOut, RayOut, LinOut

def readSK_vertexOld(fname, bigfile=False, fast_rational=False, nformat='%.14f', compresslevel=3):
    """
    Reads in Stevens vertex analysis file and returns, even more optimized for large datasets than the original.

     - a list of vertex vectors
     - a list of ray vectors
     - the basis of the lineality space as a list of vectors

    all vectors in terms of the column space of N

    """
    import gzip
    if fast_rational:
        pass
    else:
        assert _HAVE_SYMPY_, 'Install Sympy for rational IO support'

    assert os.path.exists(fname), 'Uhm exqueeze me ...'

    print('\n**********\nreadSK_vertex options are:\n')
    print('bigfile: {}'.format(bigfile))
    print('fast_rational: %s'.format(fast_rational))
    print('nformat: {}'.format(nformat % 0.12345678901234567890))
    print('**********\n')

    SK_vert_file = file(fname, 'r')
    VertOut = []
    LinOut = []
    RayOut = []
    if bigfile:
        VertTmp = gzip.open('_vtx_.tmp.gz','wb', compresslevel=compresslevel)
        LinTmp = gzip.open('_lin_.tmp.gz','wb', compresslevel=compresslevel)
        RayTmp = gzip.open('_ray_.tmp.gz','wb', compresslevel=compresslevel)
    GOvert = False
    GOray = False
    GOlin = False
    lcntr = 0
    lcntrtmp = 0
    TZero = time.time()
    for l in SK_vert_file:
        lcntr += 1
        if lcntr == 1000:
            print('Processing vertex: {} ({} min)'.format(lcntr + lcntrtmp, round((time.time()-TZero)/60.0,1)))
            lcntrtmp += lcntr
            lcntr = 0
        if '* Lineality basis ' in l:
            GOvert = False
            GOray = False
            GOlin = True
        if '* Rays ' in l:
            GOvert = False
            GOray = True
            GOlin = False
        if '* Vertices ' in l:
            GOvert = True
            GOray = False
            GOlin = False
        if l[:2] != '* ':
            L = l.split()
            rowL = []
            for c in L:
                rnum = None
                c = c.strip()
                if c == '0':
                    rnum = '0'
                    rowL.append(0.0)
                else:
                    if not fast_rational:
                        rnum = sympy.Rational('%s' % c)
                        rowL.append(rnum.evalf())
                        ##  print c, rnum.evalf()
                    else:
                        rnum = c.split('/')
                        if len(rnum) == 1:
                            rowL.append(float(rnum[0]))
                            ##  print c, float(rnum[0])
                        else:
                            rowL.append(float(rnum[0])/float(rnum[1]))
                            ##  print c, float(rnum[0].strip())/float(rnum[1].strip())
                del rnum
            del L
            rowL = tuple(rowL)
            if GOlin:
                if bigfile:
                    rowEnd = len(rowL)
                    cntr = 0
                    for e in rowL:
                        cntr += 1
                        if e == 0.0:
                            LinTmp.write('0.0')
                        else:
                            LinTmp.write(nformat % e)
                        if cntr == rowEnd:
                            LinTmp.write('\n')
                        else:
                            LinTmp.write(',')
                else:
                    LinOut.append(rowL)
            elif GOray:
                if bigfile:
                    rowEnd = len(rowL)
                    cntr = 0
                    for e in rowL:
                        cntr += 1
                        if e == 0.0:
                            RayTmp.write('0.0')
                        else:
                            RayTmp.write(nformat % e)
                        if cntr == rowEnd:
                            RayTmp.write('\n')
                        else:
                            RayTmp.write(',')
                else:
                    RayOut.append(rowL)
            elif GOvert:
                if bigfile:
                    rowEnd = len(rowL)
                    cntr = 0
                    for e in rowL:
                        cntr += 1
                        if e == 0.0:
                            VertTmp.write('0.0')
                        else:
                            VertTmp.write(nformat % e)
                        if cntr == rowEnd:
                            VertTmp.write('\n')
                        else:
                            VertTmp.write(',')
                else:
                    VertOut.append(rowL)
            del rowL
    print('\nProcessed {} vertices.\n'.format(lcntr + lcntrtmp))
    SK_vert_file.close()
    if bigfile:
        VertTmp.close()
        RayTmp.close()
        LinTmp.close()
        try:
            VertTmp = gzip.open('_vtx_.tmp.gz','rb')
            LinTmp = gzip.open('_lin_.tmp.gz','rb')
            RayTmp = gzip.open('_ray_.tmp.gz','rb')
            return VertTmp, RayTmp, LinTmp
        except Exception as ex:
            print(ex)
            print('\nReturning file names:')
            return '_vtx_.tmp.gz', '_lin_.tmp.gz', '_ray_.tmp.gz'
    else:
        print('Lineality basis: {}'.format(len(LinOut)))
        print('Number of rays: {}'.format(len(RayOut)))
        print('Number of vertices: {}'.format(len(VertOut)))
        return VertOut, RayOut, LinOut

def readSK_vertex(fname, bigfile=True, fast_rational=False, nformat='%.14f', compression=None, hdf5file=None):
    """
    Reads in Stevens vertex analysis file:

     - *fname* the input filename (.all file that results from Stevens pipeline)
     - *bigfile* [default=True] this option is now always true and is left in for backwards compatability
     - *fast_rational* [default=False] by default off and uses SymPy for rational-->float conversion, when on uses float decomposition with a slight (2th decimal) decrease in accuracy
     - *nformat* [default='%.14f'] the number format used in output files
     - *compression* [default=None] compression to be used in hdf5 files can be one of [None, 'lzf', 'gz?', 'szip']
     - *hdf5file* [default=None] if None then generic filename '_vtx_.tmp.hdf5' is uses otherwise <hdf5file>.hdf5

    and returns an hdf5 *filename* of the results with a single group named **data** which countains datasets

     - vertices
     - rays
     - lin

    where all vectors are in terms of the column space of N.

    """

    bigfile=True

    if not fast_rational:
        assert _HAVE_SYMPY_, 'Install Sympy for rational IO support'
    if bigfile:
        assert _HAVE_HD5_, 'Install h5py for large dataset support'

    assert os.path.exists(fname), 'Uhm exqueeze me ...'

    print('\n**********\nreadSK_vertex options are:\n')
    print('bigfile: {}'.format(bigfile))
    print('fast_rational: {}'.format(fast_rational))
    print('nformat: {}'.format(nformat % 0.12345678901234567890))
    print('**********\n')

    SK_vert_file = file(fname, 'r')
    VCNTR = 0
    LCNTR = 0
    RCNTR = 0
    CCNTR = 0
    for l in SK_vert_file:
        if '* Vertices (' in l:
            VCNTR = long(l.replace('* Vertices (','').replace('vectors):','').strip())
            break
    SK_vert_file.seek(0)
    for l in SK_vert_file:
        if '* Lineality basis (' in l:
            LCNTR = long(l.replace('* Lineality basis (','').replace('vectors):','').strip())
            break
    SK_vert_file.seek(0)
    for l in SK_vert_file:
        if '* Rays (' in l:
            RCNTR = long(l.replace('* Rays (','').replace('vectors):','').strip())
            break
    SK_vert_file.seek(0)
    for l in SK_vert_file:
        if l[:2] != '* ':
            L = l.split()
            CCNTR = len(L)
            del L
            break
    SK_vert_file.seek(0)
    print(RCNTR, LCNTR, VCNTR,  CCNTR)
    VertOut = []
    LinOut = []
    RayOut = []
    outFileName = '_vtx_.tmp.hdf5'
    if bigfile:
        if hdf5file != None:
            outFileName = hdf5file+'.'+str(compression)+'.hdf5'
        HD5out = h5py.File(outFileName,'w')
        if 'data' in HD5out:
            del HD5out['data']
        Dgrp = HD5out.create_group('data')
        if VCNTR > 0:
            VertTmp = HD5out['data'].create_dataset('vertices', (VCNTR, CCNTR), dtype=numpy.double, compression=compression)
        if LCNTR > 0:
            LinTmp = HD5out['data'].create_dataset('lin', (LCNTR, CCNTR), dtype=numpy.double, compression=compression)
        if RCNTR > 0:
            RayTmp = HD5out['data'].create_dataset('rays', (RCNTR, CCNTR), dtype=numpy.double, compression=compression)
    GOvert = False
    GOray = False
    GOlin = False
    lcntr = 0
    lcntrtmp = 0
    TZero = time.time()
    vert_count = 0
    lin_count = 0
    ray_count = 0
    ##  data_row = numpy.zeros((1,CCNTR), 'd')
    print('\nStarting vertex mapping at {}\n'.format(time.strftime('%H:%M:%S')))
    for l in SK_vert_file:
        lcntr += 1
        if lcntr == 1000:
            print('Processing vertex: {}: {} percent @ {} minutes ({} estimated)'.format(lcntr + lcntrtmp, (float(lcntr + lcntrtmp)/float(VCNTR)*100.0), round((time.time()-TZero)/60.0,1), float(VCNTR)/float(lcntr + lcntrtmp)*round((time.time()-TZero)/60.0,1)))
            lcntrtmp += lcntr
            lcntr = 0
        if '* Lineality basis ' in l and LCNTR > 0:
            GOvert = False
            GOray = False
            GOlin = True
        if '* Rays ' in l and RCNTR > 0:
            GOvert = False
            GOray = True
            GOlin = False
        if '* Vertices ' in l and VCNTR > 0:
            GOvert = True
            GOray = False
            GOlin = False
        if l[:2] != '* ':
            ##  data_row[0,:] = 0.0
            L = l.split()
            if GOlin:
                LinTmp[lin_count] = 0.0
                for c in xrange(CCNTR):
                    rnum = None
                    val = L[c].strip()
                    if val == '0':
                        ##  LinTmp[lin_count,c] = 0.0
                        pass
                    else:
                        if not fast_rational:
                            rnum = sympy.Rational('%s' % val)
                            LinTmp[lin_count,c] = rnum.evalf()
                        else:
                            rnum = val.split('/')
                            if len(rnum) == 1:
                                LinTmp[lin_count,c] = float(rnum[0])
                            else:
                                LinTmp[lin_count,c] = float(rnum[0])/float(rnum[1])
                lin_count += 1
            elif GOray:
                RayTmp[ray_count] = 0.0
                for c in xrange(CCNTR):
                    rnum = None
                    val = L[c].strip()
                    if val == '0':
                        ##  RayTmp[ray_count,c] = 0.0
                        pass
                    else:
                        if not fast_rational:
                            rnum = sympy.Rational('%s' % val)
                            RayTmp[ray_count,c] = rnum.evalf()
                        else:
                            rnum = val.split('/')
                            if len(rnum) == 1:
                                RayTmp[ray_count,c] = float(rnum[0])
                            else:
                                RayTmp[ray_count,c] = float(rnum[0])/float(rnum[1])
                ray_count += 1
            elif GOvert:
                VertTmp[vert_count] = 0.0
                for c in xrange(CCNTR):
                    ##  print 'lin_count',lin_count,c
                    ##  print 'ray_count',ray_count,c
                    ##  print 'vert_count',vert_count,c
                    rnum = None
                    try:
                        val = L[c].strip()
                    except Exception as ex:
                        print(ex)
                        print('\nError reading data: your file may be corrupt\n')
                    if val == '0':
                        ##  VertTmp[vert_count,c] = 0.0
                        pass
                    else:
                        if fast_rational:
                            rnum = val.split('/')
                            if len(rnum) == 1:
                                VertTmp[vert_count,c] = float(rnum[0])
                            else:
                                VertTmp[vert_count,c] = float(rnum[0])/float(rnum[1])
                        else:
                            rnum = sympy.Rational('%s' % val)
                            VertTmp[vert_count,c] = rnum.evalf()

                vert_count += 1

    print('\nProcessed {} vectors.\n'.format(lcntr + lcntrtmp))
    SK_vert_file.close()
    if bigfile:
        HD5out.close()
        ##  HD5out = h5py.File(outFileName,'r')
        return outFileName
    else:
        print('Lineality basis: {}'.format(len(LinOut)))
        print('Number of rays: {}'.format(len(RayOut)))
        print('Number of vertices: {}'.format(len(VertOut)))
        return VertOut, RayOut, LinOut

def readExcel97Model(xlname,  write_sbml=True, sbml_level=3, return_dictionaries=False):
    """
    Reads a model encoded as an Excel97 workbook and returns it as a CBMPy model object and SBML file. Note the workbook must be formatted
    exactly like those produced by cbm.writeModelToExcel97(). Note that reactions have to be defined in **both** the *reaction*
    and *network_react* sheets to be included in the model.

     - *xlpath* the filename of the Excel workbook
     - *return_model* [default=True] construct and return the CBMPy model
     - *write_sbml* [default=True] write the SBML file to fname
     - *return_dictionaries* [default=False] return the dictionaries constructed when reading the Excel file (in place of the model)
     - *sbml_level* [default=3] write the SBML file as either SBML L2 FBA or SBML L3 FBC file.


    """

    if not _HAVE_XLRD_:
        print('\nERROR: Cannot read Excel file, XLRD package not available (http://pypi.python.org/pypi/xlrd)')
        return

    assert os.path.exists(xlname), '\nERROR: File "{}" does not exist'.format(xlpath)

    def c2s(c):
        """
        Utility function converting a XLRD cell to a string
        """
        #return str(c.value.strip())
        return str(c.value)

    MSGLog = csio.StringIO()
    def logMsg(msg):
        """
        Message logging utility

        """
        print(msg)
        MSGLog.write('{}\n'.format(msg))


    wb = xlrd.open_workbook(xlname)

    for s in wb.sheets():
        logMsg('Sheet: {}'.format(s.name))

    sI = wb.sheet_by_name('info')
    sM = wb.sheet_by_name('metabolites')
    sR = wb.sheet_by_name('reactions')
    sNR = wb.sheet_by_name('network_react')
    sNM = wb.sheet_by_name('network_metab')
    sS = wb.sheet_by_name('solution')
    sMI = wb.sheet_by_name('miriam')
    try:
        sComp = wb.sheet_by_name('compartments')
    except xlrd.XLRDError:
        sComp = None


    # INFO
    dInfo = {}
    for r_ in range(3):
        r = sI.row(r_)
        if r_ == 2:
            obj = {}
            obj['osense'] = c2s(r.pop(1))
            coeff = True
            cout = []
            for c in range(1, len(r)):
                if r[c].value != '' and c2s(r[c]) != '':
                    if coeff:
                        cout.append([r[c].value])
                        coeff = False
                    else:
                        cout[-1].append(c2s(r[c]))
                        coeff = True
            obj['objflux'] = cout
            dInfo['objective'] = obj
            del cout, obj
        else:
            dInfo[c2s(r[0])] = c2s(r[1])

    # metabolites
    mcolNames = [c2s(c) for c in sM.row(0)]

    dCompartments = {}

    dMet = {}
    for r_ in range(1, sM.nrows):
        r = sM.row(r_)
        met = {}
        annot = {}
        for c_ in range(1, len(mcolNames)):
            if c_ < 6:
                if mcolNames[c_] == 'charge':
                    met[mcolNames[c_]] = r[c_].value
                elif mcolNames[c_] == 'fixed':
                    met[mcolNames[c_]] = bool(r[c_].value)
                else:
                    met[mcolNames[c_]] = c2s(r[c_])
                if mcolNames[c_] == 'compartment':
                    dCompartments[c2s(r[c_])] = {'id' : c2s(r[c_])}
            else:
                annot[mcolNames[c_]] = c2s(r[c_])
        met['annot'] = annot
        met['id'] = c2s(r[0])
        dMet[c2s(r[0])] = met

    # reactions
    rcolNames = [c2s(c) for c in sR.row(0)]

    dReact = {}
    for r_ in range(1, sR.nrows):
        r = sR.row(r_)
        reac = {}
        annot = {}
        for c_ in range(1, len(rcolNames)):
            if c_ < 6:
                if rcolNames[c_] == 'lowerbound' or rcolNames[c_] == 'upperbound':
                    reac[rcolNames[c_]] = r[c_].value
                elif rcolNames[c_] == 'reversible':
                    rev = None
                    if r[c_].ctype == 1:
                        rev = c2s(r[c_])
                        if rev in ['True', 'Yes', True, 'yes', 'true']:
                            rev = True
                        elif rev in ['False', 'No', False, 'no', 'True']:
                            rev = False
                    else:
                        rev = bool(r[c_].value)
                    reac[rcolNames[c_]] = rev
                else:
                    reac[rcolNames[c_]] = c2s(r[c_])
                if rcolNames[c_] == 'compartment':
                    dCompartments[c2s(r[c_])] = {'id' : c2s(r[c_])}
            else:
                annot[rcolNames[c_]] = c2s(r[c_])
        reac['annot'] = annot
        reac['id'] = c2s(r[0])
        dReact[c2s(r[0])] = reac

    # reaction network
    nDict = {}
    reagentList = []
    for n_ in range(0, sNR.nrows, 3):
        coeff = True
        sub = []
        s = sNR.row(n_+1)
        for c in range(1, len(s)):
            if s[c].value != '' and c2s(s[c]) != '':
                if coeff:
                    sub.append([s[c].value])
                    coeff = False
                else:
                    sub[-1].append(c2s(s[c]))
                    coeff = True
                    if c2s(s[c]) not in reagentList:
                        reagentList.append(c2s(s[c]))

        coeff = True
        prod = []
        p = sNR.row(n_+2)
        for c in range(1, len(p)):
            if p[c].value != '' and c2s(p[c]) != '':
                if coeff:
                    prod.append([p[c].value])
                    coeff = False
                else:
                    prod[-1].append(c2s(p[c]))
                    coeff = True
                    if c2s(p[c]) not in reagentList:
                        reagentList.append(c2s(p[c]))
        nDict[c2s(sNR.row(n_)[0])] = {'subs' : sub,
                                      'prod' : prod
                                      }

    # compartments
    if sComp != None:

        # metabolites
        ccolNames = [c2s(c) for c in sComp.row(0)]

        for m_ in range(1, sComp.nrows):
            r = sComp.row(m_)
            comp = {}
            annot = {}
            for c_ in range(1, len(ccolNames)):
                if c_ < 6:
                    if ccolNames[c_] == 'name':
                        comp[ccolNames[c_]] = c2s(r[c_])
                    elif ccolNames[c_] == 'size':
                        comp[ccolNames[c_]] = r[c_].value
                    elif ccolNames[c_] == 'dimensions':
                        comp[ccolNames[c_]] = r[c_].value
                    elif ccolNames[c_] == '# species' or ccolNames[c_] == '# reactions':
                        pass
                else:
                    annot[ccolNames[c_]] = c2s(r[c_])
            comp['annot'] = annot
            comp['id'] = c2s(r[0])
            dCompartments[c2s(r[0])] = comp

    # MIRIAM
    mirNames = [c2s(c) for c in sMI.row(0)]
    dMir = {}
    for m_ in range(1, sMI.nrows):
        r = sMI.row(m_)
        qualD = {}
        for c_ in range(1, len(mirNames)):
            if c2s(r[c_]) != '':
                if mirNames[c_] in qualD:
                    qualD[mirNames[c_]].append(c2s(r[c_]))
                else:
                    qualD[mirNames[c_]] = [c2s(r[c_])]
        dMir[c2s(r[0])] = qualD

    # clean up
    if '' in dCompartments:
        dCompartments.pop('')
    if '' in nDict:
        nDict.pop('')
    if '' in dReact:
        dReact.pop('')
    if '' in dMet:
        dMet.pop('')


    # cross check reactions
    dReactError = {}
    nDictError = {}
    for r_ in tuple(dReact):
        if r_ not in nDict:
            logMsg('ERROR: Reaction "{}" listed but not defined ... removing.'.format(r_))
            dReactError[r_] = dReact.pop(r_)
    for r_ in tuple(nDict):
        if r_ not in dReact:
            logMsg('ERROR: Reaction "{}" defined but not listed ... removing.'.format(r_))
            nDictError[r_] = nDict.pop(r_)

    # cross check reagents versus species
    nDictReagentError = {}
    reagentError = []

    for s_ in dMet:
        if s_ not in reagentList:
            logMsg('WARNING: Metabolite "{}" listed but not used.'.format(s_))
            #dMetError[s_] = dMet.pop(s_)

    for s_ in (reagentList):
        if s_ not in dMet:
            logMsg('ERROR: Reagent "{}" used but not defined.'.format(s_))
            reagentError.append(reagentList.pop(reagentList.index(s_)))

    # remove reactions that contain unknown reagents (undefined species)
    if len(reagentError) > 0:
        for e_ in tuple(nDict):
            reag = [a[1] for a in nDict[e_]['subs']] + [a[1] for a in nDict[e_]['prod']]
            for r_ in reag:
                if r_ in reagentError:
                    logMsg('ERROR: Reaction "{}" contains an unknown reagent "{}" ... removing'.format(e_, r_))
                    nDictReagentError[e_] = { 'reaction' : dReact.pop(e_),
                                              'stoich' : nDict.pop(e_)
                                              }
    if return_dictionaries:
        return (dInfo, dReact, dMet, nDict, dCompartments, dMir)

    # assemble CBMPy model

    # construct model object
    cmod = CBModel.Model(dInfo['id'])
    cmod.name = dInfo['name']
    cmod.description = 'Model created from Excel spreadsheet: {}'.format(xlname)
    cmod.sourcefile = xlname
    cmod.setCreatedDate()

    # add compartments
    for c_ in dCompartments:
        CC = dCompartments[c_]
        if 'name' in CC:
            name = CC['name']
        else:
            name = CC['id']
        if 'size' in CC:
            size = CC['size']
        else:
            size = 1.0
        if 'dimensions' in CC:
            dimensions = CC['dimensions']
        else:
            dimensions = 3
        C = CBModel.Compartment(CC['id'], name, size, dimensions)
        if 'annot' in CC:
            for a_ in CC['annot']:
                C.setAnnotation(a_, CC['annot'][a_])
        cmod.addCompartment(C)
    del c_, CC, name, size, dimensions

    # add metabolites
    for m_ in dMet:
        M = dMet[m_]
        S = CBModel.Species(M['id'], boundary=M['fixed'], name=M['name'], value=0.0,\
                            compartment=M['compartment'], charge=M['charge'], chemFormula=M['chemformula'])
        for a_ in M['annot']:
            S.setAnnotation(a_, M['annot'][a_])
        cmod.addSpecies(S)
    del m_, M, S

    # add reactions
    for r_ in nDict:
        Rn = nDict[r_]
        Ri = dReact[r_]
        R = CBModel.Reaction(Ri['id'], Ri['name'], reversible=Ri['reversible'])
        for s_ in Rn['subs']:
            R.createReagent(s_[1], -s_[0])
        for p_ in Rn['prod']:
            R.createReagent(p_[1], p_[0])
        for a_ in Ri['annot']:
            R.setAnnotation(a_, Ri['annot'][a_])
        cmod.addReaction(R, create_default_bounds=False)
    # cross reference reagents (shortcut)
    for s_ in cmod.species:
        tmp = s_.isReagentOf()
    del r_, s_, p_, R, Rn, Ri, tmp

    # add flux bounds
    for r_ in dReact:
        lb = dReact[r_]['lowerbound']
        ub = dReact[r_]['upperbound']
        try:
            lb = float(lb)
        except:
            if dReact[r_]['reversible']:
                lb = -numpy.inf
            else:
                lb = 0.0
            logMsg('Undefined lower bound for reaction "{}" setting to {}'.format(r_, lb))
        try:
            ub = float(ub)
        except:
            ub = numpy.inf
            logMsg('Undefined upper bound for reaction "{}" setting to {}'.format(r_, ub))
        cmod.createReactionLowerBound(r_, lb)
        cmod.createReactionUpperBound(r_, ub)
    del r_, lb, ub

    # add objective function
    O = CBModel.Objective('obj1', dInfo['objective']['osense'])
    O.createFluxObjectives(dInfo['objective']['objflux'])
    cmod.addObjective(O, True)
    del O


    # markup with MIRIAM
    for x_ in cmod.species+cmod.reactions+cmod.compartments:
        xid = x_.getId()
        if xid in dMir:
            for q_ in dMir[xid]:
                for i_ in dMir[xid][q_]:
                    #ent, mid = i_.rsplit('/',1)
                    #x_.addMIRIAMannotation(q_, ent, mid)
                    #if q_ == 'is':
                        #qual = 'isA'
                    #else:
                        #qual = q_
                    if x_.miriam == None:
                        x_.miriam = CBModel.MIRIAMannotation()
                    x_.miriam.addIDorgURI(q_, i_)

    # ok lets play
    #try:
        #cmod.createGeneAssociationsFromAnnotations()
        #geneerrors = cmod.testGeneProteinAssociations()
        #logMsg('Successfully created gene associations from annotations.')
    #except Exception as ex:
        #logMsg(ex)
    #try:
        #cbm.analyzeModel(cmod)
        #logMsg('Successfully optimized model.')
    #except Exception as ex:
        #logMsg(ex)
    #try:
        #cbm.writeModelToExcel97(cmod, xlname.replace('.xls','')+'.new')
        #logMsg('Successfully wrote model to new Excel spreadsheet: "{}"'.format(xlname.replace('.xls','')+'.DEBUG.new.xls'))
    #except Exception as ex:
        #logMsg(ex)
    if write_sbml:
        try:
            if sbml_level == 3:
                CBXML.sbml_writeSBML3FBC(cmod, xlname+'.l3.xml', gpr_from_annot=False)
                logMsg('Successfully wrote model SBML3FBC file: "{}"'.format(xlname+'.l3.xml'))
            else:
                CBXML.sbml_writeSBML2FBA(cmod, xlname+'.l3.xml')
                logMsg('Successfully wrote model SBML2FBA file: "{}"'.format(xlname+'.l3.xml'))
        except Exception as ex:
            logMsg(ex)

    # write log to file
    """
    F = file(os.path.join(cDir,'{}.log'.format(xlname)),'w')
    MSGLog.seek(0)
    F.write(MSGLog.read())
    F.close()
    """
    print('\n*****\nExcel97 Read Log\n*****\n')
    MSGLog.seek(0)
    print(MSGLog.read())
    print('')
    MSGLog.close()

    return cmod

