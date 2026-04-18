"""
CBMPy: CBTools module
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
Last edit: $Author: bgoli $ ($Id: CBTools.py 710 2020-04-27 14:22:34Z bgoli $)

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
import re
import pprint
import gzip
import zipfile

try:
    import pickle
except ImportError:
    import cPickle as pickle
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import numpy

from . import CBModel
from .CBCommon import (
    HAVE_PYPARSING,
    checkChemFormula,
    pp_chemicalFormula,
    extractGeneIdsFromString,
)
from .CBCommon import processSpeciesChargeChemFormulaAnnot, pyparsing

_PPR_ = pprint.PrettyPrinter()


from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


def createTempFileName():
    """Create a temporary filename based on the current timestamp.

    Returns a string representing a temporary filename.
    """
    return str(time.time()).split('.')[0]


# TODO comrpess


def storeObj(obj, filename, compress=False):
    """Store a Python object as a serialized binary object.

    Stores a Python *obj* as a serialized binary object in *filename*.dat.

    - *obj* a Python object
    - *filename* the base filename (will append .dat if not present)
    - *compress* [default=False] use gzip compression (not implemented)

    Returns None (writes object to disk).
    """
    if filename[-4:] != '.dat':
        filename = filename + '.dat'
    F = open(filename, 'wb')
    pickle.dump(obj, F, protocol=2)
    print('Object serialized as {}'.format(filename))
    F.close()


def loadObj(filename):
    """Load a serialized Python pickle from disk.

    Loads a serialized Python pickle from *filename*.dat and returns the Python object.

    - *filename* the base filename (will append .dat if not present)

    Returns the loaded Python object.
    """
    if filename[-4:] != '.dat':
        filename = filename + '.dat'
    assert os.path.exists(filename), '\nFile "{}" does not exist'.format(filename)
    F = open(filename, 'rb')
    obj = pickle.load(F)
    F.close()
    return obj


def deSerialize(s):
    """Deserialize a serialized object from a string.

    Deserializes a serialized Python object that was previously stored in a string.

    - *s* a string containing a serialized object

    Returns the deserialized Python object.
    """
    return pickle.loads(s)


def deSerializeFromDisk(filename):
    """Deserialize a Python pickle from a file on disk.

    Loads a serialized Python pickle from *filename* and returns the Python object.

    - *filename* the path to the serialized file

    Returns the deserialized Python object.
    """
    assert os.path.exists(filename), '\nFile "{}" does not exist'.format(filename)
    F = open(filename, 'rb')
    obj = pickle.load(F)
    F.close()
    return obj


def addStoichToFBAModel(fm):
    """Build stoichiometry for an FBA model.

    - *fm* an FBA model object

    Note: This method has been refactored into the model class. Use cmod.buildStoichMatrix() instead.

    Returns None (calls fm.buildStoichMatrix()).
    """
    fm.buildStoichMatrix()


def addSinkReaction(fbam, species, lb=0.0, ub=1000.0):
    """Add a sink reaction that consumes a model species.

    Adds a sink reaction that consumes a model *species* (X --> sink).

    - *fbam* an FBA model object
    - *species* a valid species name
    - *lb* lower flux bound [default=0.0]
    - *ub* upper flux bound [default=1000.0]

    Returns None (adds the sink reaction to the model).
    """
    assert species in fbam.getSpeciesIds(), '\n%s is not a valid species' % species
    if lb < 0.0:
        reversible = True
    else:
        reversible = False
    Rname = species + '_sink'
    R = CBModel.Reaction(
        Rname, name='%s sink reaction' % species, reversible=reversible
    )
    Su = CBModel.Reagent(Rname + species, species, -1.0)
    R.addReagent(Su)
    R.is_exchange = True

    clb = CBModel.FluxBound(Rname + '_lb', Rname, 'greaterEqual', lb)
    cub = CBModel.FluxBound(Rname + '_ub', Rname, 'lessEqual', ub)

    fbam.addReaction(R, create_default_bounds=False)
    fbam.addFluxBound(clb)
    fbam.addFluxBound(cub)

    print(
        '\n***\nCreated new reaction {} with bounds ({} : {})\n***\n'.format(
            Rname, lb, ub
        )
    )


# TODO: check this
def addSourceReaction(fbam, species, lb=0.0, ub=1000.0):
    """Add a source reaction that produces a model species.

    Adds a source reaction that produces a model *species* (source --> X).

    - *fbam* an FBA model object
    - *species* a valid species name
    - *lb* lower flux bound [default=0.0]
    - *ub* upper flux bound [default=1000.0]

    Note: Reversibility is determined by the lower bound (default 0 = irreversible).
    If negative, the reaction is reversible.

    Returns None (adds the source reaction to the model).
    """
    assert species in fbam.getSpeciesIds(), '\n%s is not a valid species' % species
    if lb < 0.0:
        reversible = True
    else:
        reversible = False
    Rname = species + '_src'
    R = CBModel.Reaction(
        Rname, name='%s source reaction' % species, reversible=reversible
    )
    Su = CBModel.Reagent(Rname + species, species, 1.0)
    R.addReagent(Su)
    R.is_exchange = True

    clb = CBModel.FluxBound(Rname + '_lb', Rname, 'greaterEqual', lb)
    cub = CBModel.FluxBound(Rname + '_ub', Rname, 'lessEqual', ub)

    fbam.addReaction(R, create_default_bounds=False)
    fbam.addFluxBound(clb)
    fbam.addFluxBound(cub)

    print(
        '\n***\nCreated new reaction {} with bounds ({} : {})\n***\n'.format(
            Rname, lb, ub
        )
    )


def findDeadEndMetabolites(fbam):
    """Find dead-end metabolites (those connected to a single reaction).

    Finds metabolites that are only involved in a single reaction (orphaned metabolites).

    - *fbam* an FBA model instance

    Returns a list of tuples (metabolite_id, reaction_id) for dead-end metabolites.
    """
    fbam.buildStoichMatrix()
    orphaned_list = []
    for rr in range(fbam.N.array.shape[0]):
        if (fbam.N.array[rr, :] != 0.0).sum() == 1:
            if __DEBUG__:
                print(fbam.N.array[rr, :])
            if __DEBUG__:
                print(fbam.N.row[rr])
            for c in range(fbam.N.array.shape[1]):
                if fbam.N.array[rr, c] != 0.0:
                    orphaned_list.append((fbam.N.row[rr], fbam.N.col[c]))
    return orphaned_list


def findDeadEndReactions(fbam):
    """Find dead-end reactions (those with a single substrate/product).

    Finds reactions that have only a single metabolite (orphaned reactions).

    - *fbam* an FBA model instance

    Returns a list of tuples (metabolite_id, reaction_id) for dead-end reactions.
    """
    fbam.buildStoichMatrix()
    orphaned_list = []
    for cc in range(fbam.N.array.shape[1]):
        if (fbam.N.array[:, cc] != 0.0).sum() == 1:
            if __DEBUG__:
                print(fbam.N.array[:, cc])
            if __DEBUG__:
                print(fbam.N.col[cc])
            for r in range(fbam.N.array.shape[0]):
                if fbam.N.array[r, cc] != 0.0:
                    orphaned_list.append((fbam.N.row[r], fbam.N.col[cc]))
    return orphaned_list


def setSpeciesPropertiesFromAnnotations(
    fbam, overwriteCharge=False, overwriteChemFormula=False
):
    """Set model species properties from annotations.

    With default options, this will only replace missing data. For ChemicalFormula
    this is easy to detect, however charge may have an "unknown value" of 0.
    Setting the optional values to True will replace any existing value with
    any valid annotation.

    - *fbam* an FBA model instance
    - *overwriteChemFormula* [default=False] whether to overwrite existing chemical formula
    - *overwriteCharge* [default=False] whether to overwrite existing charge

    Returns None (modifies species properties in place).
    """
    for s_ in fbam.species:
        try:
            processSpeciesChargeChemFormulaAnnot(
                s_,
                getFromName=False,
                overwriteCharge=overwriteCharge,
                overwriteChemFormula=overwriteChemFormula,
            )
        except Exception:
            print(
                'processSpeciesChargeChemFormulaAnnot failed for species with id: {}'.format(
                    s_.getId()
                )
            )


def fixReversibility(fbam, auto_correct=False):
    """Fix reaction reversibility by setting appropriate flux bounds.

    Set flux bound lower bounds from reactions' reversibility information.

    - *fbam* an FBA model instance
    - *auto_correct* [default=False] if True, automatically sets lower bound to zero if required.
      Otherwise, prints a warning if inconsistency is found.

    Returns None (modifies model in place).
    """
    for c in fbam.flux_bounds:
        R = c.reaction
        # print R
        O = c.operation
        # print O
        V = c.value
        # print V
        R_obj = fbam.reactions[fbam.getReactionIds().index(c.reaction)]
        RE = R_obj.reversible
        # print RE
        if O in ['greater', 'greaterEqual']:
            if not RE and float(V) < 0.0:
                print(
                    'Warning {} is not reversible and lower bound is {}.'.format(R, V)
                )
                if auto_correct:
                    print('Resetting {} lower bound ({}) to zero'.format(R, V))
                    c.value = 0.0
                else:
                    print(
                        'Reaction ({}) reversible={} inconsistent with fluxbound lower bound ({}) run with auto_correct=True to reset lower bound.'.format(
                            R, RE, V
                        )
                    )
                    time.sleep(1)


def splitReversibleReactions(fba, selected_reactions=None):
    """Split reversible reactions into forward and reverse reactions.

    Splits a set of reversible reactions into irreversible forward and reverse
    reactions (e.g., R1: A = B becomes R1f: A -> B and R1r: B -> A).

    - *fba* an instantiated CBMPy model object
    - *selected_reactions* [default=None] if provided, only split these reactions;
      otherwise splits all reversible reactions

    Returns a clone of the original model with reversible reactions split.
    """
    if selected_reactions is None:
        selected_reactions = []
    M = fba.clone()
    if len(selected_reactions) == 0:
        selected_reactions = M.getReversibleReactionIds()
    for r_ in M.getReversibleReactionIds():
        if r_ in selected_reactions:
            splitSingleReversibleReaction(M, r_)
        else:
            pass
    return M


def splitSingleReversibleReaction(fba, rid, fwd_id=None, rev_id=None):
    """Split a single reversible reaction into two irreversible reactions.

    Splits a reversible reaction (A = B) into two irreversible reactions:
    - Forward: A -> B
    - Reverse: B -> A

    The original reversible reaction and its bounds are deleted from the model.

    - *fba* an instantiated CBMPy model object
    - *rid* a valid reaction ID
    - *fwd_id* [default=None] the new forward reaction ID (defaults to rid + '_fwd')
    - *rev_id* [default=None] the new reverse reaction ID (defaults to rid + '_rev')

    Returns a tuple (original_reaction, lower_bound, upper_bound, equality_bound).
    """
    R = fba.getReaction(rid)

    assert R is not None
    print('Reversible reaction splitter is processing: {}'.format(rid))
    Rf = R.clone()
    Rb = R.clone()
    Rf.setAnnotation('cbmpy_split_fwd', rid)
    Rb.setAnnotation('cbmpy_split_rev', rid)
    RB = fba.getReactionBounds(rid)
    LB = UB = EB = None
    if RB[1] is not None and RB[2] is not None:
        assert RB[1] <= RB[2], 'ERROR: For reaction splitting ({}): LowerBound ({}) must be LessEqual to UpperBound ({})'.format(
            rid, round(RB[1], 6), round(RB[2], 6)
        )
    if RB[1] is not None:
        LB = fba.getFluxBoundByReactionID(rid, 'lower')
    if RB[2] is not None:
        UB = fba.getFluxBoundByReactionID(rid, 'upper')
    if RB[3] is not None:
        EB = fba.getFluxBoundByReactionID(rid, 'equality')
    fba.deleteReactionAndBounds(rid)
    if fwd_id is None:
        fwd_id = Rf.getId() + '_fwd'
    Rf.setPid(fwd_id)
    if rev_id is None:
        rev_id = Rb.getId() + '_rev'
    Rb.setPid(rev_id)
    Rf.reversible = False
    Rb.reversible = False
    for rr_ in Rf.reagents:
        rr_.setPid(rr_.getId() + '_fwd')
    for rr_ in Rb.reagents:
        rr_.setCoefficient(-1.0 * rr_.getCoefficient())
        rr_.setPid(rr_.getId() + '_rev')
    fba.addReaction(Rf, create_default_bounds=False)
    fba.addReaction(Rb, create_default_bounds=False)

    if EB is not None:
        fba.createReactionLowerBound(Rf.getId(), EB.getValue())
        fba.createReactionUpperBound(Rf.getId(), EB.getValue())
        fba.createReactionLowerBound(Rb.getId(), EB.getValue())
        fba.createReactionUpperBound(Rb.getId(), EB.getValue())
    elif LB is not None and UB is not None:
        if LB.getValue() <= 0.0 and UB.getValue() >= 0.0:
            fba.createReactionLowerBound(Rf.getId(), 0.0)
            fba.createReactionUpperBound(Rf.getId(), UB.getValue())
            fba.createReactionLowerBound(Rb.getId(), 0.0)
            fba.createReactionUpperBound(Rb.getId(), abs(LB.getValue()))
        elif LB.getValue() > 0.0 and UB.getValue() > 0.0:
            fba.createReactionLowerBound(Rf.getId(), LB.getValue())
            fba.createReactionUpperBound(Rf.getId(), UB.getValue())
            fba.createReactionLowerBound(Rb.getId(), 0.0)
            fba.createReactionUpperBound(Rb.getId(), 0.0)
        if LB.getValue() < 0.0 and UB.getValue() < 0.0:
            fba.createReactionLowerBound(Rf.getId(), 0.0)
            fba.createReactionUpperBound(Rf.getId(), 0.0)
            fba.createReactionLowerBound(Rb.getId(), abs(UB.getValue()))
            fba.createReactionUpperBound(Rb.getId(), abs(LB.getValue()))
    elif LB is not None and UB is None:
        if LB.getValue() > 0:
            fba.createReactionLowerBound(Rf.getId(), LB.getValue())
            fba.createReactionUpperBound(Rf.getId(), float('inf'))
            fba.createReactionLowerBound(Rb.getId(), 0.0)
            fba.createReactionUpperBound(Rb.getId(), 0.0)
        else:
            fba.createReactionLowerBound(Rf.getId(), 0.0)
            fba.createReactionUpperBound(Rf.getId(), float('inf'))
            fba.createReactionLowerBound(Rb.getId(), 0.0)
            fba.createReactionUpperBound(Rb.getId(), abs(LB.getValue()))
    elif LB is None and UB is not None:
        if UB.getValue() >= 0:
            fba.createReactionLowerBound(Rf.getId(), 0.0)
            fba.createReactionUpperBound(Rf.getId(), UB.getValue())
            fba.createReactionLowerBound(Rb.getId(), 0.0)
            fba.createReactionUpperBound(Rb.getId(), float('inf'))
        else:
            fba.createReactionLowerBound(Rf.getId(), 0.0)
            fba.createReactionUpperBound(Rf.getId(), 0.0)
            fba.createReactionLowerBound(Rb.getId(), abs(UB.getValue()))
            fba.createReactionUpperBound(Rb.getId(), float('inf'))
    else:
        fba.createReactionLowerBound(Rf.getId(), 0.0)
        fba.createReactionUpperBound(Rf.getId(), float('inf'))
        fba.createReactionLowerBound(Rb.getId(), 0.0)
        fba.createReactionUpperBound(Rb.getId(), float('inf'))

    return (R, LB, UB, EB)


def exportLabelledArray(arr, fname, names=None, sep=',', fmt='%f'):
    """Export a 2D array to a file with optional row labels.

    Write a 2D array-like object to file with optional row labels and column separator.

    - *arr* a 2D array-like object
    - *names* [default=None] the list of row names
    - *fname* the output filename
    - *sep* [default=','] the column separator character
    - *fmt* [default='%f'] the output number format

    Returns None (writes array to file).
    """
    if names is not None:
        assert arr.shape[0] == len(names), '\n ...  rows must equal number of names!'
    F = open(fname, 'w')
    cntr = 0
    for r in range(arr.shape[0]):
        if names is not None:
            F.write(('%s' + sep) % names[r])
        for c in range(arr.shape[1]):
            if c < arr.shape[1] - 1:
                F.write((fmt + sep) % arr[r, c])
            else:
                F.write((fmt + '\n') % arr[r, c])
        cntr += 1
        if cntr >= 250:
            F.flush()
            cntr = 1
    F.write('\n')
    F.flush()
    F.close()
    print('exported to {}'.format(fname))


def exportLabelledArrayWithHeader(
    arr, fname, names=None, header=None, sep=',', fmt='%f'
):
    """Export an array to file with row names and header.

    Export a 2D array to file with row labels, column header, and optional
    append mode.

    - *arr* a 2D array-like object
    - *fname* the output filename
    - *names* [default=None] the list of row names
    - *header* [default=None] the list of column names
    - *sep* [default=','] the column separator character
    - *fmt* [default='%f'] the output number format

    Returns None (writes array to file).
    """
    if names is not None:
        assert arr.shape[0] == len(names), '\n ...  rows must equal number of names!'
    if header is not None:
        assert arr.shape[1] == len(header), '\n ...  cols must equal number of header names!'

    F = open(fname, 'w')
    cntr = 0
    if header is not None:
        if names is not None:
            hstr = ' ' + sep
        else:
            hstr = ''
        for h in header:
            hstr += str(h) + sep
        hstr = hstr[:-1] + '\n'
        F.write(hstr)
        del hstr
    for r in range(arr.shape[0]):
        if names is not None:
            F.write(('%s' + sep) % names[r])
        for c in range(arr.shape[1]):
            if c < arr.shape[1] - 1:
                F.write((fmt + sep) % arr[r, c])
            else:
                F.write((fmt + '\n') % arr[r, c])
        cntr += 1
        if cntr >= 250:
            F.flush()
            cntr = 1
    F.write('\n')
    F.flush()
    F.close()
    print('exported to {}'.format(fname))


def exportLabelledLinkedList(
    arr, fname, names=None, sep=',', fmt='%s', appendlist=False
):
    """Export a 2D linked list to file with optional row labels.

    Write a 2D linked list [[...],[...],[...]] and optionally a list of row labels to file.

    - *arr* a 2D linked list (list of lists)
    - *fname* the output filename
    - *names* [default=None] the list of row names
    - *sep* [default=','] the column separator character
    - *fmt* [default='%s'] the output format string
    - *appendlist* [default=False] if True, append to file; otherwise create a new file

    Returns None (writes array to file).
    """
    if names != None:
        assert len(arr) == len(names), '\n ...  rows must equal number of names!'
    if not appendlist:
        F = open(fname, 'w')
    else:
        F = open(fname, 'a')
    cntr = 0
    for r in range(len(arr)):
        if names != None:
            F.write(('%s' + sep) % names[r])
        col_l = len(arr[0])
        for c in range(col_l):
            if c < col_l - 1:
                if arr[r][c] == 0.0:
                    F.write('0.0' + sep)
                else:
                    try:
                        F.write((fmt + sep) % arr[r][c])
                    except UnicodeEncodeError:
                        F.write((fmt + sep) % 'uError')
            else:
                if arr[r][c] == 0.0:
                    F.write('0.0\n')
                else:
                    try:
                        F.write((fmt + '\n') % arr[r][c])
                    except UnicodeEncodeError:
                        F.write((fmt + '\n') % 'uError')
        cntr += 1
        if cntr >= 250:
            F.flush()
            cntr = 1
    # F.write('\n')
    F.flush()
    F.close()
    del arr
    if not appendlist:
        print('exported to {}'.format(fname))


def exportLabelledArrayWithHeader2CSV(arr, fname, names=None, header=None):
    """Export an array with row names and header to a CSV file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .csv)
    - *names* [default=None] the list of row names
    - *header* [default=None] the list of column names

    Returns None (creates a CSV file).
    """
    fname += '.csv'
    exportLabelledArrayWithHeader(arr, fname, names, header, sep=',', fmt='%f')


def exportLabelledArray2CSV(arr, fname, names=None):
    """Export an array with row names to a CSV file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .csv)
    - *names* [default=None] the list of row names

    Returns None (creates a CSV file).
    """
    fname += '.csv'
    exportLabelledArray(arr, fname, names, sep=',', fmt='%f')


def exportArray2CSV(arr, fname):
    """Export an array to a CSV file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .csv)

    Returns None (creates a CSV file).
    """
    fname += '.csv'
    exportLabelledArray(arr, fname, None, sep=',', fmt='%f')


def exportLabelledArrayWithHeader2TXT(arr, fname, names=None, header=None):
    """Export an array with row names and header to a text file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .txt)
    - *names* the list of row names
    - *header* the list of column names

    Returns None (creates a text file).
    """
    fname += '.txt'
    exportLabelledArrayWithHeader(arr, fname, names, header, sep='\t', fmt='%f')


def exportLabelledArray2TXT(arr, fname, names=None):
    """Export an array with row names to a text file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .txt)
    - *names* [default=None] the list of row names

    Returns None (creates a text file).
    """
    fname += '.txt'
    exportLabelledArray(arr, fname, names, sep='\t', fmt='%f')


def exportArray2TXT(arr, fname):
    """Export an array to a text file.

    - *arr* a 2D array-like object
    - *fname* the output filename (will append .txt)

    Returns None (creates a text file).
    """
    fname += '.txt'
    exportLabelledArray(arr, fname, None, sep='\t', fmt='%f')


def stringReplace(fbamod, old, new, target):
    """Replace strings in model components.

    Replace occurrences of *old* with *new* in model identifiers. Target can be:

    - 'species'
    - 'reactions'
    - 'constraints'
    - 'objectives'
    - 'all' (applies to all of the above)

    Returns the modified model (fbamod).
    """
    print('stringReplace is relatively new and UNTESTED')
    fbamod.id = fbamod.id.replace(old, new)
    if target in ['species', 'all']:
        for s in fbamod.species:
            s.id = s.id.replace(old, new)
    if target in ['reactions', 'all']:
        for s in fbamod.reactions:
            s.id = s.id.replace(old, new)
            for r in s.reagents:
                r.id = r.id.replace(old, new)
    if target in ['constraints', 'all']:
        for s in fbamod.flux_bounds:
            s.id = s.id.replace(old, new)
            s.reaction = s.reaction.replace(old, new)
    if target in ['objectives', 'all']:
        for s in fbamod.objectives:
            s.id = s.id.replace(old, new)
            for f in s.flux_objectives:
                f.id = f.id.replace(old, new)
                f.reaction = f.reaction.replace(old, new)
    return fbamod


def getBoundsDict(fbamod, substring=None):
    """Get a dictionary of reactions and their bounds.

    Returns a dictionary mapping reaction names to their lower bound, upper bound,
    and equality constraints.

    - *fbamod* an FBA model instance
    - *substring* [default=None] filter reactions by ID substring

    Returns a dictionary {reaction_name: {'lb': lower_bound, 'ub': upper_bound, 'eq': equality_bound}}.
    """
    rBdic = {}
    for r in fbamod.getReactionIds(substring=substring):
        name, lb, ub, eq = fbamod.getReactionBounds(r)
        rBdic.update({name: {'lb': lb, 'ub': ub, 'eq': eq}})
    return rBdic


def getExchBoundsDict(fbamod):
    """Get bounds dictionary for all exchange reactions.

    Returns a dictionary of all exchange reactions (as determined by the is_exchange
    attribute of Reaction) after removing non-exchange reactions.

    - *fbamod* a CBMPy model

    Returns a dictionary mapping exchange reaction names to their bounds.
    """
    rBdic = {}
    for r in fbamod.getReactionIds(substring=None):
        name, lb, ub, eq = fbamod.getReactionBounds(r)
        rBdic.update({name: {'lb': lb, 'ub': ub, 'eq': eq}})
    for r in fbamod.reactions:
        if not r.is_exchange:
            rBdic.pop(r.getId())
    return rBdic


def processBiGGchemFormula(fba):
    """Disambiguate overloaded BiGG species names (NAME_CHEMFORMULA).

    Parse BiGG-style overloaded names and separate them into species name
    and chemical formula components.

    - *fba* an FBA model instance

    Note: Modifies species.name and species.chemFormula in place.
    """
    for s in fba.species:
        tmp = s.name
        tmp2 = tmp.split('_')
        if len(tmp2) >= 2:
            CF = tmp2.pop(-1)
            NM = ''
            for se in tmp2:
                NM += '%s_' % se
            NM = NM[:-1]
        else:
            NM = s.name
            CF = ''
        if __DEBUG__:
            print(NM, CF)
        del tmp, tmp2
        if s.chemFormula in ['', None, ' '] and CF != '':
            s.chemFormula = CF.strip()
            s.name = NM.strip()


def processBiGGannotationNote(fba, annotation_key='note'):
    """Parse HTML formatted reaction information from BiGG notes field.

    This function is being deprecated in favor of processSBMLAnnotationNotes().

    - *fba* an FBA model instance
    - *annotation_key* [default='note'] the annotation key containing BiGG HTML

    Warning: This function is deprecated. Use CBTools.processSBMLAnnotationNotes() instead.

    Note: Modifies r.annotation in place.
    """
    print(
        '\nDeprecation warning:\nCBTools.processBiGGannotationNote() is being replaced with CBTools.processSBMLAnnotationNotes'
    )

    html_p = re.compile("<html:p>.*?</html:p>")
    for r in fba.reactions:
        new_ann = {}
        if annotation_key in r.annotation:
            hPs = re.findall(html_p, r.annotation.pop(annotation_key))
            if __DEBUG__:
                print(hPs)
            for p in hPs:
                ps = (
                    p.replace('<html:p>', '')
                    .replace('</html:p>', '')
                    .replace('&lt;', '<')
                    .replace('&gt;', '>')
                    .split(':', 1)
                )
                if len(ps) == 2:
                    new_ann.update({ps[0].strip(): ps[1].strip()})
            r.annotation.update(new_ann)
        if __DEBUG__:
            print(r.annotation)


def processSBMLAnnotationNotes(fba, annotation_key='note', level=3):
    """Parse HTML formatted reaction information from SBML notes field.

    Parse HTML formatted reaction information stored in the SBML notes field.
    Processes BiGG and PySCeSCBM style annotations, looking for entries indexed
    with the *annotation_key*.

    - *fba* an FBA model instance
    - *annotation_key* [default='note'] the annotation key containing HTML/XHTML
      fragment (ignored in SBML Level 3)
    - *level* [default=3] SBML level (only level >= 3 uses SBML notes)

    Note: Modifies r.annotation in place.
    """

    # if hasattr(fba, '_SBML_LEVEL_') and fba._SBML_LEVEL_ != None:
    # print('\n==================================\nINFO \"CBTools.processSBMLAnnotationNotes()\":\n')
    # print('This function is now called automatically\nduring model load and can be ignored.')
    # print('==================================\n')
    # return

    html_p = re.compile("<p>.*?</p>")
    html_span = re.compile("<span>.*?</span>")
    html_bigg_p = re.compile("<html:p>.*?</html:p>")
    for r in fba.reactions:
        if level >= 3 or annotation_key in r.annotation:
            new_ann = {}
            notes = ''
            if level >= 3:
                notes = r.getNotes()
            else:
                notes = r.annotation.pop(annotation_key)
            if '<span xmlns="http://www.w3.org/1999/xhtml">' in notes:
                hPs = re.findall(html_p, notes.replace('\n', ''))
                if __DEBUG__:
                    print(hPs)
                for p in hPs:
                    ps = re.findall(html_span, p)
                    ps = [
                        p.replace('<span>', '')
                        .replace('</span>', '')
                        .replace('&lt;', '<')
                        .replace('&gt;', '>')
                        .strip()
                        for p in ps
                    ]
                    if len(ps) == 2 and ps[0] not in r.annotation:
                        new_ann.update({ps[0]: ps[1]})
            else:
                hPs = re.findall(html_bigg_p, notes)
                if len(hPs) > 0:
                    if __DEBUG__:
                        print(hPs)
                    for p in hPs:
                        ps = (
                            p.replace('<html:p>', '')
                            .replace('</html:p>', '')
                            .replace('&lt;', '<')
                            .replace('&gt;', '>')
                            .split(':', 1)
                        )
                        if len(ps) == 2 and ps[0].strip() not in r.annotation:
                            new_ann.update({ps[0].strip(): ps[1].strip()})
                else:
                    hPs = re.findall(html_p, notes)
                    if __DEBUG__:
                        print(hPs)
                    for p in hPs:
                        ps = (
                            p.replace('<p>', '')
                            .replace('</p>', '')
                            .replace('&lt;', '<')
                            .replace('&gt;', '>')
                            .split(':', 1)
                        )
                        if len(ps) == 2 and ps[0].strip() not in r.annotation:
                            new_ann.update({ps[0].strip(): ps[1].strip()})
            r.annotation.update(new_ann)
            if __DEBUG__:
                print(r.annotation)
    for s in fba.species:
        if level >= 3 or annotation_key in s.annotation:
            notes = ''
            if level >= 3:
                notes = s.getNotes()
            else:
                notes = s.annotation.pop(annotation_key)
            new_ann = {}
            if '<span xmlns="http://www.w3.org/1999/xhtml">' in notes:
                hPs = re.findall(html_p, notes.replace('\n', ''))
                if __DEBUG__:
                    print(hPs)
                for p in hPs:
                    ps = re.findall(html_span, p)
                    ps = [
                        p.replace('<span>', '')
                        .replace('</span>', '')
                        .replace('&lt;', '<')
                        .replace('&gt;', '>')
                        .strip()
                        for p in ps
                    ]
                    if len(ps) == 2 and ps[0].strip() not in s.annotation:
                        new_ann.update({ps[0]: ps[1]})
            else:
                hPs = re.findall(html_bigg_p, notes)
                if len(hPs) > 0:
                    if __DEBUG__:
                        print(hPs)
                    for p in hPs:
                        ps = (
                            p.replace('<html:p>', '')
                            .replace('</html:p>', '')
                            .replace('&lt;', '<')
                            .replace('&gt;', '>')
                            .split(':', 1)
                        )
                        if len(ps) == 2 and ps[0].strip() not in s.annotation:
                            new_ann.update({ps[0].strip(): ps[1].strip()})
                else:
                    hPs = re.findall(html_p, notes)
                    if __DEBUG__:
                        print(hPs)
                    for p in hPs:
                        ps = (
                            p.replace('<p>', '')
                            .replace('</p>', '')
                            .replace('&lt;', '<')
                            .replace('&gt;', '>')
                            .split(':', 1)
                        )
                        if len(ps) == 2 and ps[0].strip() not in s.annotation:
                            new_ann.update({ps[0].strip(): ps[1].strip()})
            s.annotation.update(new_ann)

        if 'chemFormula' in s.annotation and (
            s.chemFormula is None or s.chemFormula == ''
        ):
            s.chemFormula = s.annotation.pop('chemFormula')
            if __DEBUG__:
                print(s.annotation)
        elif 'FORMULA' in s.annotation and (
            s.chemFormula is None or s.chemFormula == ''
        ):
            s.chemFormula = s.annotation.pop('FORMULA')

        if s.chemFormula != '' and not checkChemFormula(s.chemFormula):
            s.chemFormula = ''

        if (
            (s.charge is None or s.charge == '' or s.charge == 0)
            and 'charge' in s.annotation
            and s.annotation['charge'] != ''
        ):
            chrg = s.annotation.pop('charge')
            try:
                s.charge = int(chrg)
            except ValueError:
                s.charge = None
                print(
                    'Invalid charge: {} defined for species {}'.format(chrg, s.getId())
                )
            if __DEBUG__:
                print(s.annotation)
        elif (
            (s.charge is None or s.charge == '' or s.charge == 0)
            and 'CHARGE' in s.annotation
            and s.annotation['CHARGE'] != ''
        ):
            chrg = s.annotation.pop('CHARGE')
            try:
                s.charge = int(chrg)
            except ValueError:
                print(
                    'Invalid charge: {} defined for species {}'.format(chrg, s.getId())
                )
                s.charge = None
            if __DEBUG__:
                print(s.annotation)


def processExchangeReactions(fba, key):
    """Extract exchange reactions from model and return dictionaries.

    Extract exchange reactions from model using *key* and return dictionaries
    of exchange and medium reactions.

    - *fba* an FBA model instance
    - *key* [default=None] substring to filter reactions by name

    Returns a tuple:
    - a dictionary of all exchange reactions (without medium reactions)
    - a dictionary of medium exchange reactions (those with negative lower bound)
    """
    # extract all exchange bounds
    if key is None:
        fexDic = getExchBoundsDict(fba)
    else:
        fexDic = getBoundsDict(fba, substring=key)

    # extract the medium (exchange fluxes that allow uptake)
    MediumAll = []
    Medium = []
    for r in fexDic:
        if fexDic[r]['lb'] < 0.0:
            MediumAll.append((r, fexDic[r]['lb'], fexDic[r]['ub']))
            Medium.append(r)
        if __DEBUG__:
            print(r, fexDic[r])
    # remove medium from bounds dictionary and place in medium dict
    mediumDic = {}
    for m in Medium:
        mediumDic.update({m: fexDic.pop(m)})
    if __DEBUG__:
        print('\nMedium')
        for m in MediumAll:
            print(m)
        print('mediumDic')
        print(mediumDic)
        print('\nr in fexDic')
        for r in mediumDic:
            print(r, r in fexDic)
    return fexDic, mediumDic


def generateInputScanReports(
    fba, exDict, mediumDict, optimal_growth_rates, wDir, tag=''
):
    """Generate scan reports for optimal growth rates.

    Generates CSV reports for optimal growth rates, exchange reactions,
    and medium composition.

    - *fba* an FBA model instance
    - *exDict* dictionary of exchange reaction bounds
    - *mediumDict* dictionary of medium reaction bounds
    - *optimal_growth_rates* dictionary mapping reaction IDs to optimal rates
    - *wDir* working directory for output files
    - *tag* [default=''] tag appended to model name

    Returns None (writes three CSV files: .medium.csv, .exchange.csv, .optima.csv).
    """
    modName = fba.sourcefile
    modName += tag
    rnames = fba.getReactionNames()
    rid = fba.getReactionIds()
    F = open(os.path.join(wDir, '%s.medium.csv' % modName), 'w')
    F.write('J, lb, ub\n')
    for r in mediumDict:
        RN = rnames[rid.index(r)]
        F.write('%s, %s, %s, %s\n' % (r, mediumDict[r]['lb'], mediumDict[r]['ub'], RN))
    F.write('\n')
    F.close()

    F = open(os.path.join(wDir, '%s.exchange.csv' % modName), 'w')
    F.write('J, lb, ub\n')
    for r in exDict:
        RN = rnames[rid.index(r)]
        F.write('%s, %s, %s, %s\n' % (r, exDict[r]['lb'], exDict[r]['ub'], RN))
    F.write('\n')
    F.close()

    F = open(os.path.join(wDir, '%s.optima.csv' % modName), 'w')
    F.write('J, lb, ub, optimum, "name"\n')
    if __DEBUG__:
        print(rnames)
        print(rid)
    for r in optimal_growth_rates:
        RN = rnames[rid.index(r)]
        F.write(
            '%s, %s, %s, %s, "%s"\n'
            % (r, exDict[r]['lb'], exDict[r]['ub'], optimal_growth_rates[r], RN)
        )
    F.write('\n')
    F.close()


def getAllReactionsAssociatedWithGene(
    fba, gene, gene_annotation_key='GENE ASSOCIATION'
):
    """Get all reactions associated with a gene via annotations.

    - *fba* an FBA model instance
    - *gene* a gene ID to search for
    - *gene_annotation_key* [default='GENE ASSOCIATION'] the annotation key containing gene associations

    Returns a list of reaction IDs associated with the given gene.
    """
    out = []
    for r in fba.reactions:
        GA = None
        if gene_annotation_key in r.annotation:
            GA = gene_annotation_key
        elif 'GENE ASSOCIATION' in r.annotation:
            GA = 'GENE ASSOCIATION'
        elif 'GENE_ASSOCIATION' in r.annotation:
            GA = 'GENE_ASSOCIATION'
        if GA is not None:
            if gene in r.annotation[GA]:
                out.append(r.getId())
    return out


def scanForReactionDuplicates(f, ignore_coefficients=False):
    """Scan for duplicate reactions with matching stoichiometry.

    This method uses a brute force approach to finding reactions with matching
    stoichiometry.

    - *f* an FBA model instance
    - *ignore_coefficients* [default=False] if True, ignore coefficients when comparing

    Returns a list of duplicate reaction pairs, where each pair contains:
    [reaction_id_1, reaction_id_2, refspecies, reaction_name_1, reaction_name_2]
    """
    duplicates = []
    for r in f.reactions:
        Rref = r.getSpeciesIds()
        Rref.sort()
        refspecies = ''
        for s in Rref:
            refspecies += '%s:' % s
        refspecies = refspecies[:-1]
        for r2 in f.reactions:
            Rtest = r2.getSpeciesIds()
            Rtest.sort()
            if Rref == Rtest and r.id != r2.id:
                if not ignore_coefficients:
                    go = True
                    for rgid in Rtest:
                        if float(r.getReagentWithSpeciesRef(rgid).coefficient) != float(
                            r2.getReagentWithSpeciesRef(rgid).coefficient
                        ):
                            go = False
                            break
                    if go:
                        dup = [r.id, r2.id]
                        dup.sort()
                        dup = dup + [refspecies, f.getReaction(dup[0]).getName(), f.getReaction(dup[1]).getName()]
                        if dup not in duplicates:
                            duplicates.append(dup)
                else:
                    dup = [r.id, r2.id]
                    dup.sort()
                    dup = dup + [refspecies, f.getReaction(dup[0]).getName(), f.getReaction(dup[1]).getName()]
                    if dup not in duplicates:
                        duplicates.append(dup)
    for d in duplicates:
        print(d)
    print('\nFound %s pairs of duplicate reactions' % len(duplicates))
    return duplicates


def countedPause(Tsec):
    """Print a countdown timer.

    - *Tsec* number of seconds to count down

    Prints a countdown timer, then returns.
    """
    print('\nPausing ...   ',)
    for t in range(Tsec, -1, -1):
        print('\b\b\b{}'.format(t), end=" ")
        time.sleep(1)
    print('\b\b{}'.format('done.'))


def addGenesFromAnnotations(fba, annotation_key='GENE ASSOCIATION', gene_pattern=None):
    """Deprecated method for adding genes from annotations.

    Add genes to the model using the definitions stored in the annotation key.

    - *fba* an FBA model object
    - *annotation_key* [default='GENE ASSOCIATION'] the annotation dictionary key that holds gene associations
    - *gene_pattern* [default=None] deprecated, no longer needed

    Warning: This method is deprecated. Use cmod.createGeneAssociationsFromAnnotations() instead.

    Note: Modifies the model in place.
    """

    print(
        '\nWARNING: CBTools.addGenesFromAnnotations IS DEPRECATED. USE cmod.createGeneAssociationsFromAnnotations() instead.\n'
    )
    fba.createGeneAssociationsFromAnnotations(
        annotation_key=annotation_key, replace_existing=True
    )


def getModelGenesPerReaction(
    fba, gene_pattern=None, gene_annotation_key='GENE ASSOCIATION'
):
    """Get the genes per reaction for a model.

    Parse BiGG-style gene annotation strings to extract gene IDs per reaction.

    - *fba* an FBA model instance
    - *gene_pattern* [default=None] deprecated, use model methods instead
    - *gene_annotation_key* [default='GENE ASSOCIATION'] the annotation key

    Note: It is advisable to use the model methods directly rather than this function.

    Returns a dictionary mapping reaction IDs to lists of gene IDs (or None if no genes).
    """

    react_gene = {}

    for r in fba.reactions:
        GA = None
        if gene_annotation_key in r.annotation:
            GA = gene_annotation_key
        elif 'GENE ASSOCIATION' in r.annotation:
            GA = 'GENE ASSOCIATION'
        elif 'GENE_ASSOCIATION' in r.annotation:
            GA = 'GENE_ASSOCIATION'
        elif 'gene_association' in r.annotation:
            GA = 'gene_association'
        elif 'gene association' in r.annotation:
            GA = 'gene association'
        if GA is not None:
            genes = extractGeneIdsFromString(r.annotation[GA])
            if len(genes) == 0:
                genes = None
            react_gene.update({r.getId(): genes})
        else:
            react_gene.update({r.getId(): None})
    return react_gene


def getReactionsPerGene(react_gene):
    """Get the reactions per gene for a gene-reaction dictionary.

    Convert a dictionary of gene-to-reaction mappings into a reverse mapping.

    - *react_gene* dictionary mapping gene IDs to lists of reaction IDs

    Returns a tuple (gene_react, genes, no_gene) where:
    - *gene_react*: dictionary mapping gene IDs to lists of reaction IDs
    - *genes*: list of gene IDs with associated reactions
    - *no_gene*: list of reaction IDs with no gene associated
    """
    gene_react = {}
    no_gene = []
    for R in react_gene:
        if react_gene[R] is None:
            print('Reaction {} has no gene associated with it'.format(R))
            no_gene.append(R)
        else:
            for G in react_gene[R]:
                if G in gene_react:
                    print('Updating gene {} with reaction {}'.format(G, R))
                    gene_react[G].append(R)
                else:
                    print('Adding gene {} to gene_react'.format(G))
                    gene_react.update({G: [R]})
    genes = list(gene_react)
    return gene_react, genes, no_gene


def removeFixedSpeciesReactions(f):
    """Remove fixed-species reactions from a model.

    This function removes reactions that only have boundary species as reactants
    and products. These are typically gene associations encoded in Manchester style.

    - *f* an FBA model object

    Returns None (removes reactions in place).
    """
    c_react = []
    for rea in f.reactions:
        lsa = numpy.array(
            [f.getSpecies(r.species_ref).is_boundary for r in rea.reagents]
        )
        if lsa.all():
            c_react.append(rea.getId())
    for r in c_react:
        f.deleteReactionAndBounds(r)


def addFluxAsActiveObjective(f, reaction_id, osense, coefficient=1):
    """Add a flux as an active objective function.

    - *f* an FBA model instance
    - *reaction_id* a reaction ID string
    - *osense* objective sense (must be 'maximize' or 'minimize')
    - *coefficient* [default=1] the objective function coefficient

    Returns None (modifies model in place).
    """
    osense = osense.lower()
    if osense == 'max':
        osense = 'maximize'
    if osense == 'min':
        osense = 'minimize'
    if osense in ['maximise', 'minimise']:
        osense = osense.replace('se', 'ze')
    assert osense in ['maximize', 'minimize'], (
        "\nosense must be 'maximize' or 'minimize', not %s" % osense
    )
    assert reaction_id in [r.getId() for r in f.reactions], (
        '\n%s is not a valid reaction' % reaction_id
    )
    n_obj = CBModel.Objective(reaction_id + '_objf', osense)
    f.addObjective(n_obj, active=True)
    n_flux_obj = CBModel.FluxObjective(
        reaction_id + '_fluxobj', reaction_id, coefficient
    )
    n_obj.addFluxObjective(n_flux_obj)


def checkReactionBalanceElemental(f, Rid=None, zero_tol=1.0e-12):
    """Check if reactions are balanced using chemical formulas.

    Check if reactions are balanced (elementally and by charge) using chemical
    formulas.

    - *f* an FBA model object
    - *Rid* [default=None] specific reaction(s) to check, or None for all
    - *zero_tol* [default=1.0e-12] floating point tolerance for zero check

    Requires PyParsing to be installed.

    Returns a dictionary mapping reaction IDs to their balance status:
    {'id': reaction_id, 'charge_balanced': bool, 'element_balanced': bool,
     'elements': dict, 'charge': float, 'stuff': list}
    """
    assert HAVE_PYPARSING, '\nPyParsing needs to be installed for this method'
    if Rid is None:
        Rid = f.getReactionIds()
    elif isinstance(Rid, list):
        pass
    else:
        Rid = [Rid]
    ROUT = {}
    RCHARGE = {}
    for rid in Rid:
        R = f.getReaction(rid)
        reagents = []
        netcharge = None
        for rr in R.reagents:
            CF = f.getSpecies(rr.species_ref).chemFormula
            chrg = f.getSpecies(rr.species_ref).charge
            if CF not in [None, '']:
                try:
                    CFP = pp_chemicalFormula.parseString(CF)
                    R2 = [(r[0], int(r[1])) for r in CFP]
                    reagents.append([rr.species_ref, rr.coefficient, CF, R2])
                except pyparsing.ParseException:
                    print('Invalid Chemical Formula ({}): {}'.format(rid, CF))
                    reagents.append([rr.species_ref, rr.coefficient, CF, None])
            else:
                reagents.append([rr.species_ref, rr.coefficient, CF, None])
            if chrg not in [None, '']:
                if netcharge is None:
                    netcharge = float(chrg) * rr.coefficient
                else:
                    netcharge += float(chrg) * rr.coefficient
        RCHARGE[rid] = netcharge
        ROUT[rid] = reagents
    Rres = {}
    for R in ROUT:
        Ed = {}
        for rr in ROUT[R]:
            if rr[3] is not None:
                for s in rr[3]:
                    if s[0] in Ed:
                        Ed.update({s[0]: Ed[s[0]] + rr[1] * s[1]})
                    else:
                        Ed.update({s[0]: rr[1] * s[1]})
        if len(Ed) > 0:
            CBAL = True
            EBAL = True
        else:
            CBAL = False
            EBAL = False

        for e in Ed:
            if abs(Ed[e]) >= zero_tol:
                EBAL = False
            if RCHARGE[R] is None or abs(RCHARGE[R]) >= zero_tol:
                CBAL = False
        Rres.update(
            {
                R: {
                    'id': R,
                    'charge_balanced': CBAL,
                    'element_balanced': EBAL,
                    'elements': Ed.copy(),
                    'charge': RCHARGE[R],
                    'stuff': ROUT[R],
                }
            }
        )
        if CBAL and EBAL:
            f.getReaction(R).is_balanced = True
        else:
            f.getReaction(R).is_balanced = False
    return Rres


def scanForUnbalancedReactions(f, output='all'):
    """Scan a model for unbalanced reactions.

    Scan for unbalanced reactions and return dictionaries of balanced and
    unbalanced reactions by charge and element balance.

    - *f* an FBA model instance
    - *output* [default='all'] one of ['all', 'charge', 'element']
      - 'charge': return all charge-unbalanced reactions
      - 'element': return all element-unbalanced reactions
      - 'all': return all unbalanced reactions

    Returns a dictionary of unbalanced reactions.
    """
    bcheck = checkReactionBalanceElemental(f)
    badD = bcheck.copy()
    out = {}
    all_balanced = {}
    charge_balanced = {}
    element_balanced = {}
    for b in bcheck:
        if bcheck[b]['charge_balanced'] and bcheck[b]['element_balanced']:
            all_balanced.update({b: badD.pop(b)})
        elif bcheck[b]['charge_balanced']:
            charge_balanced.update({b: badD.pop(b)})
        elif bcheck[b]['element_balanced']:
            element_balanced.update({b: badD.pop(b)})
    if output == 'charge':
        out.update(element_balanced)
    elif output == 'element':
        out.update(charge_balanced)
    else:
        out.update(element_balanced)
        out.update(charge_balanced)
    print(len(bcheck), len(badD))
    return out


def createZipArchive(zipname, files, move=False, compression='normal'):
    """Create a zip archive containing one or more files.

    - *zipname* the fully qualified name of the zip archive to create
    - *files* either a filename or list of filenames (fully qualified)
    - *move* [default=False] if True, delete input files after creating archive
    - *compression* [default='normal'] compression mode (None for no compression)

    Creates a zip archive containing one or more files.

    Returns None (writes archive to zipname).
    """

    if compression is None:
        compression = zipfile.ZIP_STORED
    else:
        compression = zipfile.ZIP_DEFLATED

    zf = zipfile.ZipFile(zipname, mode='w', compression=compression)
    if isinstance(files, list) or isinstance(files, tuple):
        files = [files]

    for f_ in files:
        assert os.path.exists(f_), 'ERROR: file \"{}\" does not exist'.format(f_)
    for f_ in files:
        zf.write(f_, arcname=os.path.split(f_)[-1])
    zf.close()
    if move:
        for f_ in files:
            try:
                os.remove(f_)
            except Exception as ex:
                print(ex)
        print(
            '\nINFO: {} input file(s) moved to archive \"{}\".'.format(
                len(files), zipname
            )
        )
    else:
        print('\nINFO: zip-archive \"{}\" created.'.format(zipname))


def checkExchangeReactions(fba, autocorrect=True):
    """Scan all reactions for exchange reactions and correct inconsistencies.

    Scan all reactions for exchange reactions (those containing a boundary species).
    Returns a list of inconsistent reactions or correct them automatically.

    - *fba* a CBMPy model
    - *autocorrect* [default=True] if True, automatically correct the "is_exchange"
      attribute on inconsistent reactions

    Returns a list of reaction IDs that were inconsistent (and possibly corrected).
    """
    badR = []
    for r_ in fba.reactions:
        has_fixed = False
        if True in [fba.getSpecies(rr_.species_ref).is_boundary for rr_ in r_.reagents]:
            has_fixed = True
        if r_.is_exchange and not has_fixed:
            print(
                'WARNING: reaction {} is labelled as an exchange reaction but has no fixed reagents.'.format(
                    r_.getId()
                )
            )
            if autocorrect:
                print('INFO: is_exchange reaction attribute corrected')
                r_.is_exchange = has_fixed
            badR.append(r_.getId())
        elif not r_.is_exchange and has_fixed:
            print(
                'WARNING: reaction {} is not labelled as an exchange reaction but contains a fixed reagent.'.format(
                    r_.getId()
                )
            )
            if autocorrect:
                print('INFO: is_exchange reaction attribute corrected')
                r_.is_exchange = has_fixed
            badR.append(r_.getId())
    return badR


def checkIds(fba, items='all'):
    """Check model attribute IDs for validity and duplicates.

    Check the IDs of specified model attributes for legal naming and duplicates.

    - *fba* a CBMPy model instance
    - *items* [default='all'] what to check: 'all' means 'species,reactions,flux_bounds,objectives'

    Returns a dictionary of errors found for each item type.
    """

    if items == 'all':
        items = [
            a.strip() for a in 'species,reactions,flux_bounds,objectives'.split(',')
        ]
    else:
        items = [a.strip() for a in items.split(',')]

    for i_ in range(len(items) - 1, -1, -1):
        if not hasattr(fba, items[i_]):
            print(
                'ERROR: bad descriptor \"{}\" removing from input list'.format(
                    items.pop(i_)
                )
            )

    output = {}
    iddump = []
    for i_ in items:
        output[i_] = []
        ITEMS = fba.__getattribute__(i_)
        for I_ in ITEMS:
            Id = I_.getId()
            if Id in iddump:
                print('INFO: duplicate \"{}\" id: {}'.format(i_, Id))
                output[i_].append(I_)
            else:
                iddump.append(Id)
            if i_ == 'reactions':
                if 'reagents' not in output:
                    output['reagents'] = []
                for rr_ in I_.reagents:
                    rrid = rr_.getId()
                    if rrid in iddump:
                        print('INFO: duplicate \"reagent\" id: {}'.format(rrid))
                        if rr_ not in output['reagents']:
                            output['reagents'].append(rr_)
                    else:
                        iddump.append(rrid)
            if i_ == 'objectives':
                if 'fluxObjectives' not in output:
                    output['fluxObjectives'] = []
                for fo_ in I_.flux_objectives:
                    foid = fo_.getId()
                    if foid in iddump:
                        print('INFO: duplicate \"fluxObjective\" id: {}'.format(foid))
                        if fo_ not in output['fluxObjectives']:
                            output['fluxObjectives'].append(fo_)
                    else:
                        iddump.append(foid)
    if len(output) == 0:
        print(
            '\nWARNING: no valid object descriptors found, please check your function call!'
        )
    return output


def checkFluxBoundConsistency(fba):
    """Check flux bound consistency for the model.

    Check for multiply defined bounds, bounds without a reaction, inconsistent
    bounds with respect to each other, and reaction reversibility.

    - *fba* a CBMPy model instance

    Returns a dictionary of errors found:
    - 'eq+lb': equality + lower bound conflicts
    - 'eq+ub': equality + upper bound conflicts
    - 'duplicate_ids': duplicate flux bound IDs
    - 'multiple_defines': multiple bound definitions
    - 'lb>ub': lower bound > upper bound
    - 'undefined': bounds without values
    - 'rev_contradict': reversibility contradictions
    - 'no_reaction': bounds without a reaction
    """

    dupIDs = checkIds(fba, items='flux_bounds')['flux_bounds']
    if len(dupIDs) > 0:
        print('\nERROR: {} duplicate flux_bound Id\'s detected!'.format(len(dupIDs)))

    LB = {}
    UB = {}
    EB = {}
    eMB = {'lower': {}, 'upper': {}, 'equality': {}}

    noreaction = []
    for fb in fba.flux_bounds:
        raw_type = fb.is_bound
        get_type = fb.getType()
        RID = fb.getReactionId()
        if raw_type != get_type:
            print(
                'WARNING: incorrect bound type for operation: \"{}\" old \"{}\" --> \"{}\"'.format(
                    fb.operation, raw_type, get_type
                )
            )

        if get_type == 'lower':
            if RID in LB:
                print(
                    'ERROR multiple LOWER bounds defined for reaction: \"{}\"'.format(
                        RID
                    )
                )
            if RID in eMB['lower']:
                eMB['lower'][RID].append(fb)
            else:
                eMB['lower'][RID] = [fb]
            LB[RID] = fb
        if get_type == 'upper':
            if RID in UB:
                print(
                    'ERROR multiple UPPER bounds defined for reaction: \"{}\"'.format(
                        RID
                    )
                )
            if RID in eMB['upper']:
                eMB['upper'][RID].append(fb)
            else:
                eMB['upper'][RID] = [fb]
            UB[RID] = fb
        if get_type == 'equality':
            if RID in EB:
                print(
                    'ERROR multiple EQUAL bounds defined for reaction: \"{}\"'.format(
                        RID
                    )
                )
            if RID in eMB['equality']:
                eMB['equality'][RID].append(fb)
            else:
                eMB['equality'][RID] = [fb]

            EB[RID] = fb
        if fba.getReaction(RID) is None:
            noreaction.append(fb)

    for mb_ in list(eMB['lower']):
        if len(eMB['lower'][mb_]) == 1:
            eMB['lower'].pop(mb_)
    for mb_ in list(eMB['upper']):
        if len(eMB['upper'][mb_]) == 1:
            eMB['upper'].pop(mb_)
    for mb_ in list(eMB['equality']):
        if len(eMB['equality'][mb_]) == 1:
            eMB['equality'].pop(mb_)

    undefined = {'no_upper': [], 'no_lower': [], 'no_upper_lower': []}

    for r_ in fba.getReactionIds():
        LBdef = True
        UBdef = True
        if r_ not in EB:
            if r_ not in LB:
                LBdef = False
            if r_ not in UB:
                UBdef = False
        if not LBdef and not UBdef:
            print('WARNING: No bounds defined for reaction: \"{}\"'.format(r_))
            undefined['no_upper_lower'].append(r_)
        else:
            if not LBdef:
                print('WARNING: No LOWER BOUND defined for reaction: \"{}\"'.format(r_))
                undefined['no_lower'].append(r_)
            if not UBdef:
                print('WARNING: No UPPER BOUND defined for reaction: \"{}\"'.format(r_))
                undefined['no_upper'].append(r_)

    errors = {
        'eq+lb': [],
        'eq+ub': [],
        'duplicate_ids': dupIDs,
        'multiple_defines': eMB,
        'lb>ub': [],
        'undefined': undefined,
        'rev_contradict': [],
        'no_reaction': noreaction,
    }

    for k_ in EB:
        if k_ in LB:
            errors['eq+lb'].append((EB[k_], LB[k_]))
        if k_ in UB:
            errors['eq+ub'].append((EB[k_], UB[k_]))

    checked = []
    for k_ in LB:
        if k_ in UB and k_ not in checked:
            if not LB[k_].getValue() <= UB[k_].getValue():
                print(
                    'ERROR: Reaction {} has lower bound ({}) larger than upper bound ({})'.format(
                        k_, LB[k_].getValue(), UB[k_].getValue()
                    )
                )
                errors['lb>ub'].append((LB[k_], UB[k_]))
            checked.append(k_)
        assR = fba.getReaction(LB[k_].getReactionId())
        if assR != None:
            if not assR.reversible:
                if LB[k_].getValue() < 0.0:
                    print(
                        'ERROR: Reaction {} is marked as irreversible but has a negative lower bound ({})'.format(
                            assR.getId(), LB[k_].getValue()
                        )
                    )
                    errors['rev_contradict'].append(assR)
            del assR

    for k_ in UB:
        if k_ in LB and k_ not in checked:
            if not LB[k_].getValue() <= UB[k_].getValue():
                print(
                    'ERROR: Reaction {} has lower bound ({}) larger than upper bound ({})'.format(
                        k_, LB[k_].getValue(), UB[k_].getValue()
                    )
                )
                errors['lb>ub'].append((LB[k_], UB[k_]))
            checked.append(k_)

    return errors


def roundOffWithSense(val, osense='max', tol=1e-8):
    """Round a value with sense-aware rounding.

    Round a value in a way that takes into consideration the sense of the
    operation that generated it (maximize rounds up, minimize rounds down).

    - *val* the value to round
    - *osense* [default='max'] the operation sense ('maximize' or 'minimize')
    - *tol* [default=1e-8] the tolerance for roundoff

    Returns the rounded value.
    """
    if osense.lower() in ['min', 'minimize', 'minimise']:
        val = numpy.ceil(val / tol) * tol
    else:
        val = numpy.floor(val / tol) * tol
    return val


def mergeGroups(m, groups, new_id, new_name='', auto_delete=False):
    """Merge a list of groups into a new group.

    Merge a list of groups into a new group. Note: annotations are not merged.

    - *m* a model containing the source groups
    - *groups* a list of group IDs to merge
    - *new_id* the new merged group ID
    - *new_name* [default=''] the new group name (default merges old names)
    - *auto_delete* [default=False] if True, delete source groups after merge

    Returns True if merge succeeded, False otherwise.
    """
    if type(groups) is list and len(groups) > 1:
        badgid = []
        m_gids = m.getGroupIds()
        for gnew in groups:
            if gnew not in m_gids:
                badgid.append(gnew)
        if len(badgid) > 0:
            print('ERROR: groups contains invalid group ids: {}'.format(str(badgid)))
            return False
    else:
        print('ERROR: groups must be a list with more than one element.')
        return False

    if m.getGroup(new_id) is not None:
        print('ERROR: new_id {} already exists'.format(new_id))
        return False

    m.createGroup(new_id)
    gnew = m.getGroup(new_id)

    make_name = False
    if new_name == '':
        make_name = True

    for gid in groups:
        gobj = m.getGroup(gid)
        if make_name:
            new_name = '{}+{}'.format(new_name, gobj.getName())
        for gm in gobj.members:
            if gm not in gnew.members:
                gnew.addMember(gm)
        if auto_delete:
            m.deleteGroup(gid)
    if make_name:
        new_name = new_name[1:]
    gnew.setName(new_name)

    return True


def merge2Models(m1, m2, ignore=None, ignore_duplicate_ids=False):
    """Merge two models (m2 into m1) without model checking.

    Merge 2 models, this method does a raw merge of model 2 into model 1
    without any model checking. Component IDs in ignore are ignored in both
    models and the first objective of model 1 is arbitrarily set as active.
    Compartments are also merged and a new "OuterMerge" compartment is created.

    In all cases, duplicate IDs are tracked and ignored, essentially using
    the object ID encountered first (usually that of model 1). Duplicate
    checking can be disabled by setting the *ignore_duplicate_ids* flag.

    - *m1* model 1 (target)
    - *m2* model 2 (source)
    - *ignore* [default=None] list of IDs to skip during merge
    - *ignore_duplicate_ids* [default=False] disable duplicate ID checking

    Note: In development - merging genes and GPRs.

    Returns the merged model.
    """
    if ignore is None:
        ignore = []

    out = CBModel.Model(m1.getId() + m2.getId())
    out.setName(m1.getName() + m2.getName())
    out.createCompartment('OuterMerge', size=1.0, dimensions=3)

    idstore = []
    for x_ in m1.compartments + m2.compartments:
        sid = x_.getId()
        if sid not in ignore:
            if ignore_duplicate_ids or sid not in idstore:
                idstore.append(sid)
                out.addCompartment(x_.clone())

    for s_ in m1.species + m2.species:
        sid = s_.getId()
        if sid not in ignore:
            if ignore_duplicate_ids or sid not in idstore:
                idstore.append(sid)
                out.addSpecies(s_.clone())
            else:
                print('Skipping duplicate id: "{}"'.format(sid))
        else:
            print('Skipping ignored id: "{}"'.format(sid))
    for r_ in m1.reactions + m2.reactions:
        sid = r_.getId()
        if r_.getId() not in ignore:
            if ignore_duplicate_ids or sid not in idstore:
                idstore.append(sid)
                out.addReaction(r_.clone(), create_default_bounds=False)
            else:
                print('Skipping duplicate id: "{}"'.format(sid))
        else:
            print('Skipping ignored id: "{}"'.format(sid))
    for f_ in m1.flux_bounds + m2.flux_bounds:
        sid = f_.getId()
        if f_.getId() not in ignore:
            if ignore_duplicate_ids or sid not in idstore:
                idstore.append(sid)
                out.addFluxBound(f_.clone())
            else:
                print('Skipping duplicate id: "{}"'.format(sid))
        else:
            print('Skipping ignored id: "{}"'.format(sid))
    GO = True
    for o_ in m1.objectives + m2.objectives:
        sid = o_.getId()
        if o_.getId() not in ignore:
            if ignore_duplicate_ids or sid not in idstore:
                idstore.append(sid)
                if GO:
                    out.addObjective(o_.clone(), active=True)
                    GO = False
                else:
                    out.addObjective(o_.clone(), active=False)
            else:
                print('Skipping duplicate id: "{}"'.format(sid))
        else:
            print('Skipping ignored id: "{}"'.format(sid))
    print('\nAdded {} components to merged model'.format(len(idstore)))
    idstore = []
    return out
