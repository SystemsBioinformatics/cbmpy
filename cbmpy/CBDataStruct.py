"""
CBMPy: CBDataStruct module
==========================
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
Last edit: $Author: bgoli $ ($Id: CBDataStruct.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, copy, re
import numpy

import webbrowser
from . import miriamids

MIRIAM = miriamids.miriamids

#MIRIAM.update({'EC' : {'data_entry': 'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&amp;ec=$id',
                       #'example': '1.1.1.1',
                       #'name': 'Enzyme Nomenclature',
                       #'pattern': re.compile(r'^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$'),
                       #'url': 'http://identifiers.org/ec-code/'}
               #})
MIRIAM_KEYS = list(MIRIAM)
MIRIAM_KEYS.sort()
MIRIAM_KEYS = tuple(MIRIAM_KEYS)
MIRIAM_KEYSlc = tuple([a.lower() for a in MIRIAM_KEYS])


from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

cTime = time.time

# work around Python version insanity
if os.sys.version_info[0] >= 3 and os.sys.version_info[1] >= 7:
    RE_PATTERN_TYPE = re.Pattern
else:
    RE_PATTERN_TYPE = re._pattern_type

class StructMatrix:
    """
    This class is specifically designed to store structural matrix information
    give it an array and row/col index permutations it can generate its own
    row/col labels given the label src.
    """

    array = None
    ridx = None
    cidx = None
    row = None
    col = None
    shape = None

    def __init__(self, array, ridx, cidx, row=None, col=None):
        """
        Instantiate with array and matching row/col index arrays, optional label arrays
        """
        self.array = array
        self.ridx = ridx
        self.cidx = cidx
        self.row = row
        self.col = col
        self.shape = array.shape

    def __call__(self):
        return self.array

    def getRowsByIdx(self, *args):
        """Return the rows referenced by index (1,3,5)"""
        return self.array.take(args, axis=0)

    def getColsByIdx(self, *args):
        """Return the columns referenced by index (1,3,5)"""
        return self.array.take(args, axis=1)

    def setRow(self, src):
        """
        Assuming that the row index array is a permutation (full/subset)
        of a source label array by supplying that source to setRow it
        maps the row labels to ridx and creates self.row (row label list)
        """
        self.row = [src[r] for r in self.ridx]

    def setCol(self, src):
        """
        Assuming that the col index array is a permutation (full/subset)
        of a source label array by supplying that src to setCol
        maps the row labels to cidx and creates self.col (col label list)
        """

        self.col = [src[c] for c in self.cidx]

    def getRowsByName(self, *args):
        """Return the rows referenced by label ('s','x','d')"""
        assert self.row != None, "\nI need row labels"
        try:
            return self.array.take([self.row.index(l) for l in args], axis=0)
        except Exception as ex:
            print(ex)
            print("\nValid row labels are: {}".format(self.row))
            return None

    def getColsByName(self, *args):
        """Return the columns referenced by label ('s','x','d')"""
        assert self.col != None, "\nI need column labels"
        try:
            return self.array.take([self.col.index(l) for l in args], axis=1)
        except Exception as ex:
            print(ex)
            print("Valid column labels are: {}".format(self.col))
            return None

    def getLabels(self, axis='all'):
        """Return the matrix labels ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row': return self.row
        elif axis == 'col': return self.col
        else: return self.row, self.col

    def getIndexes(self, axis='all'):
        """Return the matrix indexes ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row': return self.ridx
        elif axis == 'col': return self.cidx
        else: return self.ridx, self.cidx

    def getByIdx(self, row, col):
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        return self.array[row, col]

    def getByName(self, row, col):
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        return self.array[self.row.index(row), self.col.index(col)]

    def setByIdx(self, row, col, val):
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        self.array[row, col] = val

    def setByName(self, row, col, val):
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        self.array[self.row.index(row), self.col.index(col)] = val

class StructMatrixLP(StructMatrix):
    """Adds some stuff to StructMatrix that makes it LP friendly"""
    RHS = None
    operators = None
    __array_type__ = None

    def __init__(self, array, ridx, cidx, row=None, col=None, rhs=None, operators=None):
        """
        Object that holds an LP, stoichiometric constraints, operators, RHS, names etc

         - *array* stoichiometric matrix (linear constraints)
         - *ridx* a list of indexes typically range(shape[0])
         - *cidx* a list of indexes typically range(shape[1])
         - *row* a list of row names
         - *col* a list of col names
         - *rhs* [default=None] a list of rhs values defaults to 0
         - *operators* [default=None] a list of constraint senses defaults to 'E' =

        """
        StructMatrix.__init__(self, array, ridx, cidx, row=row, col=col)
        self.__array_type__ = type(array)
        if type(rhs) == type(None):
            self.RHS = numpy.zeros(array.shape[0])
        else:
            assert len(rhs) == array.shape[0], "\nRHS length mismatch"
            if type(rhs) == list:
                rhs = numpy.array(rhs)
            self.RHS = rhs
        if operators == None:
            self.operators = ['E']*len(self.RHS)
        else:
            assert len(operators) == len(self.RHS), "\nOperator length mismatch"
            self.operators = operators

    def getCopy(self, attr_str, deep=False):
        """
        Return a copy of the attribute with name attr_str. Uses the copy module `copy.copy` or `copy.deepcopy`

        - *attr_str* a string of the attribute name: 'row', 'col'
        - *deep* [default=False] try to do a deepcopy. Use with caution see copy module docstring for details

        """
        if attr_str == 'array' or attr_str == 'rhs':
            return getattr(self, attr_str).copy()
        elif deep:
            return copy.deepcopy(getattr(self, attr_str))
        else:
            return copy.copy(getattr(self, attr_str))

    def getLandRHS(self):
        RHS = self.RHS.copy()
        RHS.resize(len(RHS), 1)
        return numpy.hstack([self.array, RHS])

    def setRHSbyName(self, name, value):
        assert name in self.row, '\n%s is an invalid name' % name
        self.RHS[self.row.index(name)] = value

    def setRHSbyIdx(self, idx, val):
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        self.RHS[idx] = val

    def setOperatorbyName(self, name, value):
        assert name in self.row, '\n%s is an invalid name' % name
        value = value.strip()
        assert value in ['G','L','E'], '\n%s is not a valid operator' % value
        self.operators[self.row.index(name)] = value

    def setOperatorbyIdx(self, idx, value):
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        value = value.strip()
        assert value in ['G','L','E'], '\n%s is not a valid operator' % value
        self.operators[idx] = value

    def getRHSbyName(self, name):
        assert name in self.row, '\n%s is an invalid name' % name
        return self.RHS[self.row.index(name)]

    def getRHSbyIdx(self, idx):
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        return self.RHS[idx]

    def getOperatorbyName(self, name):
        assert name in self.row, '\n%s is an invalid name' % name
        return self.operators[self.row.index(name)]

    def getOperatorbyIdx(self, idx):
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        return self.operators[idx]

#TODO: this all needs to be redone to allow for nested qualifiers and model/biol qualifers

class MIRIAMannotation(object):
    """
    The MIRIAMannotation class MIRIAM annotations: Biological Qualifiers

    """

    MIRIAM = None
    MIDS = None
    MIDSlc = None
    QUALIFIERS = None
    # biomodels qualifiers
    isA = None
    isEncodedBy = None
    encodes = None
    hasPart = None
    hasProperty = None
    hasTaxon = None
    hasVersion = None
    isDescribedBy = None
    isHomologTo = None
    isPartOf = None
    isPropertyOf = None
    isVersionOf = None
    occursIn = None
    isDerivedFrom = None

    def __init__(self):
        self.MIRIAM = MIRIAM
        self.MIDS = MIRIAM_KEYS
        self.MIDSlc = MIRIAM_KEYSlc
        self.QUALIFIERS = ("isA","isEncodedBy","encodes","hasPart","hasProperty","hasTaxon","hasVersion","isDescribedBy",\
                           "isHomologTo","isPartOf","isPropertyOf","isVersionOf","occursIn",\
                           "isDerivedFrom")

    def addIDorgURI(self, qual, uri):
        """
        Add a URI directly into a qualifier collection:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *uri* the complete identifiers.org uri e.g. http://identifiers.org/chebi/CHEBI:58088

        """
        if qual == 'is':
            qual = 'isA'
        if hasattr(self, qual):
            Q = self.__getattribute__(qual)
            if Q == None:
                Q = []
            Q.append(uri)
            self.__setattr__(qual, Q)
        else:
            print('INFO: Invalid qualifier: \"{}\" uri NOT set'.format(qual))
        #print(self.getAllMIRIAMUris())

    def checkEntityPattern(self, entity):
        """
        For an entity key compile the pattern to a regex, if necessary.

         - *entity* a MIRIAM resource entity

        """

        if type(self.MIRIAM[entity]['pattern']) == RE_PATTERN_TYPE:
            return True
        else:
            try:
                #self.MIRIAM[entity]['pattern'] = re.compile(self.MIRIAM[entity]['pattern'])
                return True
            except:
                print('Invalid pattern (entity={}): \"{}\"'.format(self.MIRIAM[entity]['name'], self.MIRIAM[entity]['pattern']))
                return False

    def checkId(self, entity, mid):
        """
        Check that a entity id e.g. CHEBI:17158

         - *mid* the entity id e.g. CHEBI:17158

        """
        if self.checkEntityPattern(entity):
            res = re.findall(self.MIRIAM[entity]['pattern'], mid)
            if len(res) > 0:
                return True
        print('INFO: invalid entity ({}) id: \"{}\"'.format(entity, mid))
        return False

    def checkEntity(self, entity):
        """
        Check an entity entry, this is a MIRIAM resource name: "chEBI". The test is case insensitive and will correct the case
        of wrongly capitalised entities automatically. If the entity is not recognised then a list of possible candidates
        based on the first letters of the input is displayed.

        - *entity* a MIRIAM resource entity e.g. "ChEBI"

        """
        if entity.lower() not in self.MIDSlc:
            print('\n\"{}\" is not a valid entity were you looking for one of these:\n'.format(entity))
            if len(entity) == 0:
                print('Need something to work with here')
            elif len(entity) == 1:
                temp = [a for a in self.MIDS if a[0].lower() == entity[0].lower()]
            else:
                temp = [a for a in self.MIDS if a[:2].lower() == entity[:2].lower()]
            for t_ in temp:
                print('\t{}'.format(t_))
            return None
        else:
            return self.MIDS[self.MIDSlc.index(entity.lower())]

    def addMIRIAMannotation(self, qual, entity, mid):
        """
        Add a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158

        """
        if qual == 'is':
            qual = 'isA'
        if hasattr(self, qual):
            Q = self.__getattribute__(qual)
            if Q == None:
                Q = []
            E = self.checkEntity(entity)
            if E != None:
                if self.checkId(E, mid):
                    Q.append(self.MIRIAM[E]['url']+mid)
                    self.__setattr__(qual, Q)
            else:
                print('INFO: Invalid entity: \"{}\" MIRIAM entity NOT set'.format(entity))
        else:
            print('INFO: Invalid qualifier: \"{}\" MIRIAM entity NOT set'.format(qual))

    def deleteMIRIAMannotation(self, qual, entity, mid):
        """
        Deletes a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158

        """
        if qual == 'is':
            qual = 'isA'
        if hasattr(self, qual):
            Q = self.__getattribute__(qual)
            if Q != None:
                E = self.checkEntity(entity)
                if E != None:
                    if self.checkId(E, mid):
                        annot = self.MIRIAM[E]['url']+mid
                        if annot in Q:
                            idx = Q.index(annot)
                            print('Deleted annotation {}'.format(Q.pop(idx)))
                        else:
                            print('INFO: Invalid entity: \"{}\" MIRIAM entry does not exist'.format(mid))
                else:
                    print('INFO: Invalid entity: \"{}\" MIRIAM entry does not exist'.format(mid))
            else:
                print('INFO: Invalid entity: \"{}\" MIRIAM entity does not exist'.format(entity))
        else:
            print('INFO: Invalid qualifier: \"{}\" MIRIAM entity does not exist'.format(qual))

    def getMIRIAMUrisForQualifier(self, qual):
        """
        Return all list of urls associated with qualifier:

         - *qual* the qualifier e.g. "is" or "isEncoded"

        """
        if qual == 'is':
            qual = 'isA'
        if hasattr(self, qual):
            Q = self.__getattribute__(qual)
            if Q == None:
                print('INFO: No \"{}\" qualifiers defined'.format(qual))
                return ()
            else:
                return tuple(Q)
        else:
            print('INFO: Invalid qualifier: \"{}\" MIRIAM entity NOT set'.format(qual))
            return ()

    def getAndViewUrisForQualifier(self, qual):
        """
        Retrieve all url's associated with qualifier and attempt to open them all in a new browser tab

         - *qual* the qualifier e.g. "is" or "isEncoded"

        """
        Q = self.getMIRIAMUrisForQualifier(qual)
        if Q != None:
            for q_ in Q:
                self.viewURL(q_)

    def getAllMIRIAMUris(self):
        """
        Return a dictionary of qualifiers that contain ID.org URL'S
        """
        out = {}
        for q_ in self.QUALIFIERS:
            Q = self.__getattribute__(q_)
            urls = []
            if Q != None:
                for u_ in self.getMIRIAMUrisForQualifier(q_):
                    urls.append(u_)
            if q_ == 'isA':
                q_ = 'is'
            out[q_] = tuple(urls)
        return out

    def viewURL(self, url):
        """
        This will try to open the URL in a new tab of the default webbrowser

         - *url* the url

        """
        webbrowser.open_new_tab(url)

#class MIRIAMModelAnnotation(MIRIAMannotation):
    #"""
    #Derived class with qualifiers for BQmodel
    #"""
    #isDerivedFrom = None


    #def __init__(self):
        #self.MIRIAM = MIRIAM
        #self.MIDS = MIRIAM_KEYS
        #self.MIDSlc = MIRIAM_KEYSlc
        #self.QUALIFIERS = ("isA", "isDerivedFrom", "isDescribedBy")



##  class SensitivityReaction(object):
    ##  """
    ##  Pseudo reaction class holding sensitivity data
    ##  """

    ##  def __init__(self):
        ##  pass

##  class SensitivitySpecies(object):
    ##  """
    ##  Pseudo species class holding sensitivity data
    ##  """

    ##  def __init__(self):
        ##  pass

##  class SensitivityData(object):
    ##  """
    ##  A class that holds sensitivity analysis data of an LP object
    ##  """

    ##  def __init__(self, sense_obj, sense_rhs, sense_bnd):
        ##  pass

    ##  def getListOfReactions(self):
        ##  """
        ##  Returns a list of sensitivity reaction objects
        ##  """
        ##  return []

    ##  def getListOfConstraints(self):
        ##  """
        ##  Returns a list of sensitivity constraint objects
        ##  """
        ##  return []
