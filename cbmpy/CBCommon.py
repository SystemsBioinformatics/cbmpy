"""
CBMPy: CBCommon module
======================
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
Last edit: $Author: bgoli $ ($Id: CBCommon.py 710 2020-04-27 14:22:34Z bgoli $)

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
import copy
import re
import numpy
import weakref
import ast
import webbrowser

from . import miriamids
from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']
MIRIAM = miriamids.miriamids
MIRIAM_KEYS = list(MIRIAM)
MIRIAM_KEYS.sort()
MIRIAM_KEYS = tuple(MIRIAM_KEYS)
MIRIAM_KEYSlc = tuple([a.lower() for a in MIRIAM_KEYS])

# TODO fine out what this comment is about
# MIRIAM.update({'EC' : {'data_entry': 'http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&amp;ec=$id',
# 'example': '1.1.1.1',
# 'name': 'Enzyme Nomenclature',
# 'pattern': re.compile(r'^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$'),
# 'url': 'http://identifiers.org/ec-code/'}
# })

# work around Python version insanity
if os.sys.version_info[0] >= 3 and os.sys.version_info[1] >= 7:
    RE_PATTERN_TYPE = re.Pattern
else:
    RE_PATTERN_TYPE = re._pattern_type

chemElementDefs = None
HAVE_PYPARSING = False

try:
    import pyparsing
    HAVE_PYPARSING = True
    print('pyparsing import')
except ImportError:
    print('INFO: please install pyparsing (pip install pyparsing) to enable chemical balance checking!')


#if not HAVE_PYPARSING:
    #if os.sys.version_info[0] == 3 and os.sys.version_info[1] < 9:
        #try:
            #from . import pyparsing

            #HAVE_PYPARSING = True
        #except ImportError:
            #print('INFO please install pyparsing to enable chemical balance checking!')

if HAVE_PYPARSING:
    pp_caps = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    pp_lowers = pp_caps.lower()
    pp_digits = "0123456789"
    pp_element = pyparsing.Word(pp_caps, pp_lowers)
    pp_integer = pyparsing.Word(pp_digits)
    pp_elementRef = pyparsing.Group(
        pp_element + pyparsing.Optional(pp_integer, default="1")
    )
    pp_chemicalFormula = pyparsing.OneOrMore(pp_elementRef)

    pp_gene_ass = (
        pyparsing.Word(pyparsing.alphanums + '. ') | 'or' | 'OR' | 'and' | 'AND'
    )
    pp_gene_unnester = pyparsing.nestedExpr(content=pp_gene_ass)
else:
    print('Please install pyparsing! pip install pyparsing')

ptElements = (
    'H',
    'He',
    'Li',
    'Be',
    'B',
    'C',
    'N',
    'O',
    'F',
    'Ne',
    'Na',
    'Mg',
    'Al',
    'Si',
    'P',
    'S',
    'Cl',
    'Ar',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'G',
    'As',
    'Se',
    'Br',
    'Kr',
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
    'Cs',
    'Ba',
    'La',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
    'Fr',
    'Ra',
    'Ac',
    'Rf',
    'Db',
    'Sg',
    'Bh',
    'Hs',
    'Mt',
    'Ds',
    'Rg',
    '112',
    '113',
    '114',
    '115',
    '116',
    '117',
    '118',
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
    'X',
    'R',
    'Z',
)


def parseGeneAssociation(gs):
    """Parse a COBRA style gene association into a nested list.

    Parameters
    ----------
    gs : str
        A string containing a gene association.

    Returns
    -------
    list
        A nested list representation of the gene association.
    """
    # gs = '(b0810) and ( b0811 ) or ( b1234.0) and(b0809)and ( b7643 )OR(b0812 )AND( b0876)'
    gl = pp_gene_unnester.parseString('(' + gs + ')').asList()
    gl = gl[0]

    def f(l):
        for i in range(len(l)):
            try:
                l[i] = l[i].lower().strip()
            except Exception:
                f(l[i])

    f(gl)
    print(gl)
    return gl


def checkChemFormula(cf, quiet=False):
    """Check if a string conforms to a Chemical Formula.

    Please see the SBML Level 3 specification and http://wikipedia.org/wiki/Hill_system
    for more information.

    Parameters
    ----------
    cf : str
        A string that contains a formula to check.
    quiet : bool, optional
        Do not print error messages. Default is False.

    Returns
    -------
    bool
        True if valid chemical formula, False otherwise.
    """
    if not HAVE_PYPARSING:
        print('\nChemical Formula parser requires PyParsing!')
        return False
    try:
        R2 = [(r[0], int(r[1])) for r in pp_chemicalFormula.parseString(cf)]
    except Exception:
        if not quiet:
            print('WARNING: \"{}\" is not a valid chemical formula.'.format(cf))
        return False

    cf2 = ''
    for l_ in R2:
        if l_[0] not in ptElements:
            # print(R2)
            if not quiet:
                print('WARNING: \"{}\" is not a valid chemical formula.'.format(cf))
            return False
        else:
            cf2 += l_[0]
            if l_[1] != 1:
                cf2 += str(l_[1])

    if cf != cf2:
        if not quiet:
            print(R2)
            print('INFO: \"{}\" check formula.'.format(cf))
        return False
    else:
        return True


def extractGeneIdsFromString(g, return_clean_gpr=False):
    """Extract and return a list of gene names from a gene association string
    formulation.

    Parameters
    ----------
    g : str
        A COBRA style gene association string.
    return_clean_gpr : bool, optional
        In addition to the list returns the "cleaned" GPR string. Default is False.

    Returns
    -------
    list or tuple
        A list of gene names, or a tuple of (names, cleaned_gpr) if return_clean_gpr=True.
    """
    # print('\n{}'.format(g))
    g2 = g
    if 'AND' in g2 or 'and' in g2:
        g2 = g2.replace(')AND ', ') AND ').replace(')and ', ') and ')
        g2 = g2.replace(' AND(', ' AND (').replace(' and(', ' and (')
        g2 = g2.replace(')AND(', ') AND (').replace(')and(', ') and (')
    if 'OR' in g2 or 'or' in g2:
        g2 = g2.replace(')OR ', ') OR ').replace(')or ', ') or ')
        g2 = g2.replace(' OR(', ' OR (').replace(' or(', ' or (')
        g2 = g2.replace(')OR(', ') OR (').replace(')or(', ') or (')

    g2 = g2.replace(' AND ', ' and ').replace(' OR ', ' or ')
    g3 = g2.replace('(', '').replace(')', '')
    tmp = []

    if ' and ' in g3 or ' or ' in g3:
        if ' and ' and not ' or ' in g3:
            tmp = g3.split(' and ')
        elif ' or ' and not ' and ' in g3:
            tmp = g3.split(' or ')
        else:
            tmp1 = g3.split(' or ')
            tmp = []
            for s_ in tmp1:
                if ' and ' in s_:
                    tmp3 = s_.split(' and ')
                    tmp += tmp3
                else:
                    tmp += [s_]
    else:
        tmp = [g3]
    names = []
    for n_ in tmp:
        n = n_.strip()
        if n not in names:
            names.append(n)
    if names == ['']:
        names = []

    if not return_clean_gpr:
        return names
    else:
        return names, g2


def createAssociationDictFromNode(node, out, model, useweakref=True, cntr=0):
    """Convert a GPR string '((g1 and g2) or g3)' to a dictionary via a Python AST.

    Parameters
    ----------
    node : ast.AST
        A Python AST node (e.g. body).
    out : dict
        A GPR dictionary (modified in place).
    model : Fbase
        A CBMPy model instance.
    useweakref : bool, optional
        Use a weakref as the gene object or alternatively the label. Default is True.
    cntr : int, optional
        Counter for unique key generation. Default is 0.

    Returns
    -------
    None
    """
    if isinstance(node, ast.Name):
        gene = model.getGene(node.id)
        if gene is not None:
            if useweakref:
                out[node.id] = weakref.ref(gene)
            else:
                out[node.id] = gene.getLabel()
        del gene

    elif isinstance(node, ast.BinOp):
        left = node.left.id
        right = node.right.id
        gref = '{}-{}'.format(left, right)
        print(gref)
        # ref.setGeneProduct(formatSbmlId(gref))
    else:
        if isinstance(node, ast.Expr):
            children = [node.value]
        else:
            children = node.values
        for v in children:
            if isinstance(v, ast.BoolOp) and isinstance(v.op, ast.And):
                key = '_AND_{}'.format(cntr)
                out[key] = {}
                newex = out[key]
                cntr += 1
                # print('And', v)
                # walk(v, newand)
            elif isinstance(v, ast.BoolOp) and isinstance(v.op, ast.Or):
                key = '_OR_{}'.format(cntr)
                out[key] = {}
                newex = out[key]
                cntr += 1
                # print('Or', v)
                # walk(v, newor)
            else:
                # print('-->', v)
                newex = out
            createAssociationDictFromNode(
                v, newex, model, useweakref=useweakref, cntr=cntr
            )


def func_getAssociationStrFromGprDict(gprd, out, parent=''):
    """Get an old school GPR association string from a CBMPy gprDict.

    Parameters
    ----------
    gprd : dict
        The GPR dictionary.
    out : str
        The output string accumulator.
    parent : str, optional
        The string representing the current node's parent relationship, used for recursion.
        Default is ''.

    Returns
    -------
    str
        The GPR association string.
    """
    out2 = '('
    for k in gprd:
        if k.startswith('_AND_'):
            out2 += '{}{}'.format(
                func_getAssociationStrFromGprDict(gprd[k], out, ' and '), parent
            )
        elif k.startswith('_OR_'):
            out2 += '{}{}'.format(
                func_getAssociationStrFromGprDict(gprd[k], out, ' or '), parent
            )
        else:
            out2 += '{}{}'.format(k, parent)
    if out2.endswith(' and '):
        out2 = out2[:-5]
    elif out2.endswith(' or '):
        out2 = out2[:-4]
    out2 += ')'
    out = out + out2
    if ' and ' not in out and ' or ' not in out:
        out = out[1:-1]
    return out


def getGPRasDictFromString(node, out):
    """Convert a GPR string '((g1 and g2) or g3)' to a gprDict.

    Parameters
    ----------
    node : ast.AST
        A Python AST node (e.g. ``ast.parse(gprstring).body[0]``).
    out : dict
        A new dictionary that will be created in place.

    Returns
    -------
    dict
        The populated gpr dictionary.
    """

    if isinstance(node, ast.Name):
        # print('Gene: {}'.format(node.id))
        out[node.id] = node.id
    elif isinstance(node, str):
        out[node] = node
    elif isinstance(node, ast.Str):
        out[node.s] = node.s
    elif isinstance(node, ast.BinOp):
        left = node.left.id
        right = node.right.id
        gref = '{}-{}'.format(left, right)
        # print('BinOp: {}'.format(gref))
    else:
        if isinstance(node, ast.Expr):
            children = [node.value]
        elif isinstance(node, ast.BoolOp):
            children = node.values
        else:
            print('Childless node: ', type(node))
            children = []
        cntr2 = 0
        for v in children:
            if isinstance(v, ast.BoolOp) and isinstance(v.op, ast.And):
                out['_AND_{}'.format(cntr2)] = {}
                # print('And', v)
                getGPRasDictFromString(v, out['_AND_{}'.format(cntr2)])
            elif isinstance(v, ast.BoolOp) and isinstance(v.op, ast.Or):
                out['_OR_{}'.format(cntr2)] = {}
                # print('Or', v)
                getGPRasDictFromString(v, out['_OR_{}'.format(cntr2)])
            else:
                # print('-->', v)
                getGPRasDictFromString(v, out)
            cntr2 += 1
    return out


class ComboGen(object):
    """Generate sets of unique combinations.

    This class generates and manages combinations of elements.

    Attributes
    ----------
    combo : list
        List of stored combinations.
    combo_int : list
        Integer representation of combinations.
    """

    combo = None
    combo_int = None

    def __init__(self):
        """Initialize the ComboGen instance.

        Returns
        -------
        None
        """
        self.combo = []

    def addCombination(self, data):
        """Add a combination to the internal list.

        Parameters
        ----------
        data : list
            A list of elements to join into a string and add to combo.

        Returns
        -------
        None
        """
        self.combo.append(','.join(data))

    def uniqueCombinations(self, data, number, temp=None):
        """Generate unique combinations of the given data.

        Parameters
        ----------
        data : list
            A list of elements to combine.
        number : int
            The number of elements per combination.
        temp : list, optional
            Temporary list used internally. Default is None.

        Returns
        -------
        None
        """
        if temp == None:
            temp = []
        if not number:
            self.addCombination(temp)
        for i in range(len(data)):
            temp.append(data[i])
            self.uniqueCombinations(data[i + 1 :], number - 1, temp)
            temp.pop()

    def numberifyComb2Int(self):
        """Convert string combinations to integer tuples.

        Returns
        -------
        None
        """
        tmp = []
        for c in self.combo:
            tmp.append(tuple([int(x) for x in c.split(',')]))
        self.combo_int = tmp


def processSpeciesChargeChemFormulaAnnot(
    s, getFromName=False, overwriteChemFormula=False, overwriteCharge=False
):
    """Disambiguate the chemical formula and charge from species annotations.

    Parameters
    ----------
    s : object
        A species object.
    getFromName : bool, optional
        Whether to try stripping the chemical formula from the name (old COBRA style).
        Default is False.
    overwriteChemFormula : bool, optional
        Whether to overwrite existing chemical formula. Default is False.
    overwriteCharge : bool, optional
        Whether to overwrite existing charge. Default is False.

    Returns
    -------
    None
    """
    tempF = ''
    if overwriteChemFormula:
        tempF = s.chemFormula
        s.chemFormula = None
    key = None
    if s.chemFormula == None or s.chemFormula == '':
        if 'chemFormula' in s.annotation:
            key = 'chemFormula'
        elif 'FORMULA' in s.annotation:
            key = 'FORMULA'
        if key != None:
            s.chemFormula = s.annotation[key]

    if not checkChemFormula(s.chemFormula, quiet=True):
        if getFromName:
            try:
                n = s.getName()
                n, cf = n.rsplit('_', 1)
                if checkChemFormula(cf, quiet=True):
                    s.chemFormula = cf
                    s.setName(n)
                else:
                    s.chemFormula = ''
            except Exception:
                s.chemFormula = ''
        else:
            s.chemFormula = ''
        if overwriteChemFormula and s.chemFormula == '':
            s.chemFormula = tempF

    else:
        if key != None:
            s.annotation.pop(key)

    key2 = None
    tempC = None
    if overwriteCharge:
        tempC = s.charge
        s.charge = None
    if s.charge == None or s.charge == '':
        if 'charge' in s.annotation:
            key2 = 'charge'
        elif 'CHARGE' in s.annotation:
            key2 = 'CHARGE'
        if key2 != None:
            chrg = s.annotation[key2]
            try:
                s.charge = int(chrg)
            except ValueError:
                print(
                    'Invalid charge: {} defined for species {}'.format(chrg, s.getId())
                )
                s.charge = None
                key2 = None
    if overwriteCharge and s.charge == None:
        s.charge = tempC
    if key2 != None:
        s.annotation.pop(key2)


def binHash(keys, d):
    """Create a hash tuple for binning by checking keys against a dictionary.

    Parameters
    ----------
    keys : list
        A list of keys to check.
    d : dict
        A dictionary.

    Returns
    -------
    tuple
        A tuple of values from d for keys present, or True for missing keys.
    """
    out = []
    for k in keys:
        if k in d:
            out.append(d[k])
        else:
            out.append(True)
    return tuple(out)


def fixId(s, replace=None):
    """Fix a string to be a valid C style variable name.

    The first letter must be an underscore or letter, and the rest should be alphanumeric or underscore.

    Parameters
    ----------
    s : str
        The string to test.
    replace : str or None, optional
        Character to use for illegal characters, or None to remove them. Default is None.

    Returns
    -------
    str
        A fixed string conforming to C variable naming rules.
    """

    if s[0].isdigit():
        s2 = '_'
        # print('Illegal start character \"{}\": prefixing an underscore'.format(s[0]))
    else:
        s2 = ''
    for c in s:
        if c.isalnum() or c == '_':
            s2 += c
        elif replace != None:
            s2 += replace
            # print('Replacing illegal character \"{}\": {} --> {}'.format(c, s, s2))
        # else:
        # print('Skipping illegal character \"{}\": {} --> {}'.format(c, s, s2))
    return s2


def checkId(s):
    """Check if a string conforms to a C variable naming convention.

    Parameters
    ----------
    s : str
        A string.

    Returns
    -------
    bool
        True if the string is a valid C variable name, False otherwise.
    """
    cntr = 0
    for c in s:
        if cntr == 0 and c.isalpha() or c == '_':
            pass
        elif cntr > 0 and c.isalnum() or c == '_':
            pass
        else:
            print('\"{}\" is an invalid character in \"{}\"'.format(c, s))
            return False
        cntr += 1
    return True


class StructMatrix:
    """Store structural matrix information and generate row/col labels.

    This class stores structural matrix data and uses row/column index permutations
    to generate labels from a source.

    Attributes
    ----------
    array : ndarray
        The structural matrix array.
    ridx : list
        Row indexes.
    cidx : list
        Column indexes.
    row : list
        Row labels.
    col : list
        Column labels.
    shape : tuple
        Shape of the matrix.
    """

    array = None
    ridx = None
    cidx = None
    row = None
    col = None
    shape = None

    def __init__(self, array, ridx, cidx, row=None, col=None):
        """Instantiate a StructMatrix with array and matching row/col index arrays.

        Parameters
        ----------
        array : ndarray
            The structural matrix array.
        ridx : list
            Array of row indexes.
        cidx : list
            Array of column indexes.
        row : list or None, optional
            Optional array of row labels. Default is None.
        col : list or None, optional
            Optional array of column labels. Default is None.

        Returns
        -------
        None
        """
        self.array = array
        self.ridx = ridx
        self.cidx = cidx
        self.row = row
        self.col = col
        self.shape = array.shape

    def __call__(self):
        """Return the stored array.

        Returns
        -------
        ndarray
            The stored array.
        """
        return self.array

    def getRowsByIdx(self, *args):
        """Return the rows referenced by index.

        Parameters
        ----------
        *args : int
            Indices of rows to retrieve (e.g. 1, 3, 5).

        Returns
        -------
        ndarray
            A subset of rows from the array.
        """
        return self.array.take(args, axis=0)

    def getColsByIdx(self, *args):
        """Return the columns referenced by index.

        Parameters
        ----------
        *args : int
            Indices of columns to retrieve (e.g. 1, 3, 5).

        Returns
        -------
        ndarray
            A subset of columns from the array.
        """
        return self.array.take(args, axis=1)

    def setRow(self, src):
        """Map row labels to ridx using a source label array.

        Assuming that the row index array is a permutation of a source label array,
        this method maps the row labels to ridx and creates self.row (row label list).

        Parameters
        ----------
        src : list
            A source array of row labels.

        Returns
        -------
        None
        """
        self.row = [src[r] for r in self.ridx]

    def setCol(self, src):
        """Map column labels to cidx using a source label array.

        Assuming that the col index array is a permutation of a source label array,
        this method maps the row labels to cidx and creates self.col (col label list).

        Parameters
        ----------
        src : list
            A source array of column labels.

        Returns
        -------
        None
        """

        self.col = [src[c] for c in self.cidx]

    def getRowsByName(self, *args):
        """Return the rows referenced by label.

        Parameters
        ----------
        *args : str
            Row labels to retrieve (e.g. 's', 'x', 'd').

        Returns
        -------
        ndarray or None
            A subset of rows from the array, or None if row labels are not set.
        """
        assert self.row != None, "\nI need row labels"
        try:
            return self.array.take([self.row.index(l) for l in args], axis=0)
        except Exception as ex:
            print(ex)
            print("\nValid row labels are: {}".format(self.row))
            return None

    def getColsByName(self, *args):
        """Return the columns referenced by label.

        Parameters
        ----------
        *args : str
            Column labels to retrieve (e.g. 's', 'x', 'd').

        Returns
        -------
        ndarray or None
            A subset of columns from the array, or None if column labels are not set.
        """
        assert self.col != None, "\nI need column labels"
        try:
            return self.array.take([self.col.index(l) for l in args], axis=1)
        except Exception as ex:
            print(ex)
            print("Valid column labels are: {}".format(self.col))
            return None

    def getLabels(self, axis='all'):
        """Return the matrix labels.

        Parameters
        ----------
        axis : str, optional
            'row'/'col'/'all' to specify which labels to return. Default is 'all'.

        Returns
        -------
        list or tuple
            A single label list or a tuple of (row_labels, col_labels).
        """
        if axis == 'row':
            return self.row
        elif axis == 'col':
            return self.col
        else:
            return self.row, self.col

    def getIndexes(self, axis='all'):
        """Return the matrix indexes.

        Parameters
        ----------
        axis : str, optional
            'row'/'col'/'all' to specify which indexes to return. Default is 'all'.

        Returns
        -------
        list or tuple
            A single index list or a tuple of (row_indexes, col_indexes).
        """
        if axis == 'row':
            return self.ridx
        elif axis == 'col':
            return self.cidx
        else:
            return self.ridx, self.cidx

    def getByIdx(self, row, col):
        """Return the array value at the given row and column index.

        Parameters
        ----------
        row : int
            A row index.
        col : int
            A column index.

        Returns
        -------
        float
            The value at self.array[row, col].
        """
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        return self.array[row, col]

    def getByName(self, row, col):
        """Return the array value at the given row and column name.

        Parameters
        ----------
        row : str
            A row name.
        col : str
            A column name.

        Returns
        -------
        float
            The value at self.array[row, col].
        """
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        return self.array[self.row.index(row), self.col.index(col)]

    def setByIdx(self, row, col, val):
        """Set the array value at the given row and column index.

        Parameters
        ----------
        row : int
            A row index.
        col : int
            A column index.
        val : float
            The value to set.

        Returns
        -------
        None
        """
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        self.array[row, col] = val

    def setByName(self, row, col, val):
        """Set the array value at the given row and column name.

        Parameters
        ----------
        row : str
            A row name.
        col : str
            A column name.
        val : float
            The value to set.

        Returns
        -------
        None
        """
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        self.array[self.row.index(row), self.col.index(col)] = val


class StructMatrixLP(StructMatrix):
    """Adds LP-friendly attributes to StructMatrix.

    This class extends StructMatrix to add LP (linear programming) support
    including RHS values, constraint operators, and stoichiometric constraints.

    Attributes
    ----------
    RHS : ndarray
        Right-hand side values for constraints.
    operators : list
        List of constraint senses ('G', 'L', 'E').
    __array_type__ : type
        Type of the underlying array.
    """

    RHS = None
    operators = None
    __array_type__ = None

    def __init__(self, array, ridx, cidx, row=None, col=None, rhs=None, operators=None):
        """Initialize a StructMatrixLP with LP-friendly attributes.

        Object that holds an LP (linear programming) model with stoichiometric constraints,
        operators, RHS values, and names.

        Parameters
        ----------
        array : ndarray
            Stoichiometric matrix (linear constraints).
        ridx : list
            A list of indexes typically range(shape[0]).
        cidx : list
            A list of indexes typically range(shape[1]).
        row : list or None, optional
            A list of row names. Default is None.
        col : list or None, optional
            A list of col names. Default is None.
        rhs : array-like or None, optional
            A list of RHS values (defaults to 0). Default is None.
        operators : list or None, optional
            A list of constraint senses (defaults to 'E'). Default is None.

        Returns
        -------
        None
        """
        StructMatrix.__init__(self, array, ridx, cidx, row=row, col=col)
        self.__array_type__ = type(array)
        if rhs is None:
            self.RHS = numpy.zeros(array.shape[0])
        else:
            assert len(rhs) == array.shape[0], "\nRHS length mismatch"
            if type(rhs) == list:
                rhs = numpy.array(rhs)
            self.RHS = rhs
        if operators is None:
            self.operators = ['E'] * len(self.RHS)
        else:
            assert len(operators) == len(self.RHS), "\nOperator length mismatch"
            self.operators = operators

    def getCopy(self, attr_str, deep=False):
        """Return a copy of the specified attribute.

        Uses the copy module ``copy.copy`` or ``copy.deepcopy`` depending on the deep flag.

        Parameters
        ----------
        attr_str : str
            A string of the attribute name (e.g. 'row', 'col').
        deep : bool, optional
            If True, uses deepcopy; otherwise uses copy. Default is False.

        Returns
        -------
        ndarray or list
            A copy of the attribute.
        """
        if attr_str == 'array' or attr_str == 'rhs':
            return getattr(self, attr_str).copy()
        elif deep:
            return copy.deepcopy(getattr(self, attr_str))
        else:
            return copy.copy(getattr(self, attr_str))

    def getLandRHS(self):
        """Return the stoichiometric matrix augmented with RHS column.

        Returns
        -------
        ndarray
            A numpy array of [self.array, RHS].
        """
        RHS = self.RHS.copy()
        RHS.resize(len(RHS), 1)
        return numpy.hstack([self.array, RHS])

    def setRHSbyName(self, name, value):
        """Set the RHS value by row name.

        Parameters
        ----------
        name : str
            A row name.
        value : float
            The RHS value to set.

        Returns
        -------
        None
        """
        assert name in self.row, '\n%s is an invalid name' % name
        self.RHS[self.row.index(name)] = value

    def setRHSbyIdx(self, idx, val):
        """Set the RHS value by row index.

        Parameters
        ----------
        idx : int
            A row index.
        val : float
            The RHS value to set.

        Returns
        -------
        None
        """
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        self.RHS[idx] = val

    def setOperatorbyName(self, name, value):
        """Set the constraint operator by row name.

        Parameters
        ----------
        name : str
            A row name.
        value : str
            A constraint operator ('G', 'L', 'E').

        Returns
        -------
        None
        """
        assert name in self.row, '\n%s is an invalid name' % name
        value = value.strip()
        assert value in ['G', 'L', 'E'], '\n%s is not a valid operator' % value
        self.operators[self.row.index(name)] = value

    def setOperatorbyIdx(self, idx, value):
        """Set the constraint operator by row index.

        Parameters
        ----------
        idx : int
            A row index.
        value : str
            A constraint operator ('G', 'L', 'E').

        Returns
        -------
        None
        """
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        value = value.strip()
        assert value in ['G', 'L', 'E'], '\n%s is not a valid operator' % value
        self.operators[idx] = value

    def getRHSbyName(self, name):
        """Get the RHS value by row name.

        Parameters
        ----------
        name : str
            A row name.

        Returns
        -------
        float
            The RHS value for the given row.
        """
        assert name in self.row, '\n%s is an invalid name' % name
        return self.RHS[self.row.index(name)]

    def getRHSbyIdx(self, idx):
        """Get the RHS value by row index.

        Parameters
        ----------
        idx : int
            A row index.

        Returns
        -------
        float
            The RHS value for the given index.
        """
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        return self.RHS[idx]

    def getOperatorbyName(self, name):
        """Get the constraint operator by row name.

        Parameters
        ----------
        name : str
            A row name.

        Returns
        -------
        str
            The operator ('G', 'L', or 'E') for the given row.
        """
        assert name in self.row, '\n%s is an invalid name' % name
        return self.operators[self.row.index(name)]

    def getOperatorbyIdx(self, idx):
        """Get the constraint operator by row index.

        Parameters
        ----------
        idx : int
            A row index.

        Returns
        -------
        str
            The operator ('G', 'L', or 'E') for the given index.
        """
        assert idx in self.ridx, '\n%s is an invalid index' % idx
        return self.operators[idx]


# TODO: this all needs to be redone to allow for nested qualifiers and model/biol qualifers


class MIRIAMannotation(object):
    """MIRIAM annotations for Biological Qualifiers.

    This class handles MIRIAM (Mapping Identification Resource) annotations,
    which provide standardized identifiers for biological entities.

    Attributes
    ----------
    MIRIAM : dict
        MIRIAM resource identifiers.
    MIDS : tuple
        MIRIAM keys as tuple.
    MIDSlc : tuple
        Lowercase MIRIAM keys as tuple.
    QUALIFIERS : tuple
        Biomodels qualifier names.
    isA : list or None
        'is' qualifier entries.
    isEncodedBy : list or None
        'isEncodedBy' qualifier entries.
    encodes : list or None
        'encodes' qualifier entries.
    hasPart : list or None
        'hasPart' qualifier entries.
    hasProperty : list or None
        'hasProperty' qualifier entries.
    hasTaxon : list or None
        'hasTaxon' qualifier entries.
    hasVersion : list or None
        'hasVersion' qualifier entries.
    isDescribedBy : list or None
        'isDescribedBy' qualifier entries.
    isHomologTo : list or None
        'isHomologTo' qualifier entries.
    isPartOf : list or None
        'isPartOf' qualifier entries.
    isPropertyOf : list or None
        'isPropertyOf' qualifier entries.
    isVersionOf : list or None
        'isVersionOf' qualifier entries.
    occursIn : list or None
        'occursIn' qualifier entries.
    isDerivedFrom : list or None
        'isDerivedFrom' qualifier entries.
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
        """Initialize the MIRIAMannotation instance.

        Returns
        -------
        None
        """
        self.MIRIAM = MIRIAM
        self.MIDS = MIRIAM_KEYS
        self.MIDSlc = MIRIAM_KEYSlc
        self.QUALIFIERS = (
            "isA",
            "isEncodedBy",
            "encodes",
            "hasPart",
            "hasProperty",
            "hasTaxon",
            "hasVersion",
            "isDescribedBy",
            "isHomologTo",
            "isPartOf",
            "isPropertyOf",
            "isVersionOf",
            "occursIn",
            "isDerivedFrom",
        )

    def addIDorgURI(self, qual, uri):
        """Add a URI directly into a qualifier collection.

        Parameters
        ----------
        qual : str
            A Biomodels biological qualifier (e.g. 'is', 'isEncodedBy').
        uri : str
            The complete identifiers.org URI (e.g. 'http://identifiers.org/chebi/CHEBI:58088').

        Returns
        -------
        None
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
        # print(self.getAllMIRIAMUris())

    def checkEntityPattern(self, entity):
        """Compile the regex pattern for a MIRIAM resource entity.

        For an entity key, this compiles the pattern to a regex if necessary.

        Parameters
        ----------
        entity : str
            A MIRIAM resource entity.

        Returns
        -------
        bool
            True if the pattern is valid/compiled, False otherwise.
        """

        if type(self.MIRIAM[entity]['pattern']) == RE_PATTERN_TYPE:
            return True
        else:
            try:
                # self.MIRIAM[entity]['pattern'] = re.compile(self.MIRIAM[entity]['pattern'])
                return True
            except:
                print(
                    'Invalid pattern (entity={}): \"{}\"'.format(
                        self.MIRIAM[entity]['name'], self.MIRIAM[entity]['pattern']
                    )
                )
                return False

    def checkId(self, entity, mid):
        """Check that a MIRIAM entity ID is valid.

        Parameters
        ----------
        entity : str
            A MIRIAM resource entity (e.g. 'CHEBI').
        mid : str
            The entity ID (e.g. 'CHEBI:17158').

        Returns
        -------
        bool
            True if the ID matches the entity pattern, False otherwise.
        """
        if self.checkEntityPattern(entity):
            res = re.findall(self.MIRIAM[entity]['pattern'], mid)
            if len(res) > 0:
                return True
        print('INFO: invalid entity ({}) id: \"{}\"'.format(entity, mid))
        return False

    def checkEntity(self, entity):
        """Check if a MIRIAM resource entity name is valid.

        This checks a MIRIAM resource name (e.g., 'chEBI'). The test is case
        insensitive and automatically corrects wrongly capitalized entities. If
        the entity is not recognized, a list of possible candidates based on
        the first letters of the input is displayed.

        Parameters
        ----------
        entity : str
            A MIRIAM resource entity (e.g. 'ChEBI').

        Returns
        -------
        str or None
            The normalized entity name if valid, None otherwise.
        """
        if entity.lower() not in self.MIDSlc:
            print(
                '\n\"{}\" is not a valid entity were you looking for one of these:\n'.format(
                    entity
                )
            )
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
        """Add a qualified MIRIAM annotation.

        Parameters
        ----------
        qual : str
            A Biomodels biological qualifier (e.g. 'is', 'isEncodedBy').
        entity : str
            A MIRIAM resource entity (e.g. 'ChEBI').
        mid : str
            The entity ID (e.g. 'CHEBI:17158').

        Returns
        -------
        None
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
                    Q.append(self.MIRIAM[E]['url'] + mid)
                    self.__setattr__(qual, Q)
            else:
                print(
                    'INFO: Invalid entity: \"{}\" MIRIAM entity NOT set'.format(entity)
                )
        else:
            print('INFO: Invalid qualifier: \"{}\" MIRIAM entity NOT set'.format(qual))

    def deleteMIRIAMannotation(self, qual, entity, mid):
        """Delete a qualified MIRIAM annotation.

        Parameters
        ----------
        qual : str
            A Biomodels biological qualifier (e.g. 'is', 'isEncodedBy').
        entity : str
            A MIRIAM resource entity (e.g. 'ChEBI').
        mid : str
            The entity ID (e.g. 'CHEBI:17158').

        Returns
        -------
        None
        """
        if qual == 'is':
            qual = 'isA'
        if hasattr(self, qual):
            Q = self.__getattribute__(qual)
            if Q != None:
                E = self.checkEntity(entity)
                if E != None:
                    if self.checkId(E, mid):
                        annot = self.MIRIAM[E]['url'] + mid
                        if annot in Q:
                            idx = Q.index(annot)
                            print('Deleted annotation {}'.format(Q.pop(idx)))
                        else:
                            print(
                                'INFO: Invalid entity: \"{}\" MIRIAM entry does not exist'.format(
                                    mid
                                )
                            )
                else:
                    print(
                        'INFO: Invalid entity: \"{}\" MIRIAM entry does not exist'.format(
                            mid
                        )
                    )
            else:
                print(
                    'INFO: Invalid entity: \"{}\" MIRIAM entity does not exist'.format(
                        entity
                    )
                )
        else:
            print(
                'INFO: Invalid qualifier: \"{}\" MIRIAM entity does not exist'.format(
                    qual
                )
            )

    def getMIRIAMUrisForQualifier(self, qual):
        """Return all URLs associated with a qualifier.

        Parameters
        ----------
        qual : str
            The qualifier (e.g. 'is' or 'isEncoded').

        Returns
        -------
        tuple
            A tuple of URLs for the given qualifier, or () if not found.
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
        """Retrieve all URLs associated with a qualifier and open them in browser tabs.

        Parameters
        ----------
        qual : str
            The qualifier (e.g. 'is' or 'isEncoded').

        Returns
        -------
        None
        """
        Q = self.getMIRIAMUrisForQualifier(qual)
        if Q is not None:
            for q_ in Q:
                self.viewURL(q_)

    def getAllMIRIAMUris(self):
        """Return a dictionary of qualifiers that contain ID.org URLs.

        Returns
        -------
        dict
            A dictionary mapping qualifier names to tuples of their associated
            ID.org URLs.
        """
        out = {}
        for q_ in self.QUALIFIERS:
            Q = self.__getattribute__(q_)
            urls = []
            if Q is not None:
                for u_ in self.getMIRIAMUrisForQualifier(q_):
                    urls.append(u_)
            if q_ == 'isA':
                q_ = 'is'
            out[q_] = tuple(urls)
        return out

    def viewURL(self, url):
        """Open a URL in a new browser tab.

        Parameters
        ----------
        url : str
            The URL to open.

        Returns
        -------
        None
        """
        webbrowser.open_new_tab(url)
