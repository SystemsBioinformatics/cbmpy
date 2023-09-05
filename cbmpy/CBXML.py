"""
CBMPy: CBXML module
===================
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

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBXML.py 710 2020-04-27 14:22:34Z bgoli $)

"""
## gets rid of "invalid variable name" info
# pylint: disable=C0103
## gets rid of "line to long" info
# pylint: disable=C0301
## use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# all unicode literals must first be cast to string for libSBML API - brett
# from __future__ import unicode_literals

import os, time, numpy, re
import ast, shutil

try:
    from cgi import escape as ESCAPE
except:
    from html import escape as ESCAPE

try:
    from cgi import escape as ESCAPE
except:
    from html import escape as ESCAPE

# this is a hack that needs to be streamlined a bit
try:
    import cStringIO as csio
except ImportError:
    import io as csio
from getpass import getuser
import xml.etree.ElementTree as ELTree
from xml.dom.minidom import getDOMImplementation

try:
    from HTMLParser import HTMLParser
except ImportError:
    # from html.parser import HTMLParser
    from html import parser

    HTMLParser = parser.HTMLParser

from . import CBModel
from .CBCommon import (
    getGPRasDictFromString,
    processSpeciesChargeChemFormulaAnnot,
    MIRIAMannotation,
)
from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

_HAVE_SBML_ = None
try:
    import libsbml

    _HAVE_SBML_ = True
    SBMLreader = libsbml.SBMLReader()
    SBMLwriter = libsbml.SBMLWriter()

    SBML_TYPES = {
        libsbml.SBML_SPECIES: 'species',
        libsbml.SBML_REACTION: 'reaction',
        libsbml.SBML_COMPARTMENT: 'compartment',
        libsbml.SBML_PARAMETER: 'parameter',
        libsbml.SBML_INITIAL_ASSIGNMENT: 'initialassignment',
        libsbml.SBML_SPECIES_REFERENCE: 'speciesreference',
        libsbml.SBML_UNKNOWN: 'unknown',
    }
    try:
        SBML_TYPES.update(
            {
                libsbml.SBML_FBC_FLUXBOUND: 'fluxbound',
                libsbml.SBML_FBC_OBJECTIVE: 'objective',
                libsbml.SBML_FBC_FLUXOBJECTIVE: 'fluxobjective',
                libsbml.SBML_FBC_V1ASSOCIATION: 'geneassociation1',
                libsbml.SBML_FBC_GENEASSOCIATION: 'geneassociation2',
                libsbml.SBML_FBC_GENEPRODUCT: 'geneproduct',
                libsbml.SBML_FBC_GENEPRODUCTASSOCIATION: 'geneproductassociation',
                libsbml.SBML_FBC_GENEPRODUCTREF: 'geneproductref',
            }
        )
        FBC3_VARIABLE_TYPES = {libsbml.FBC_VARIABLE_TYPE_LINEAR : 'linear',
                              libsbml.FBC_VARIABLE_TYPE_QUADRATIC : 'quadratic',
                              libsbml.FBC_VARIABLE_TYPE_INVALID : 'invalid',
                              'linear' : libsbml.FBC_VARIABLE_TYPE_LINEAR,
                              'quadratic' : libsbml.FBC_VARIABLE_TYPE_QUADRATIC,
                              'invalid' : libsbml.FBC_VARIABLE_TYPE_INVALID
                              }
    except AttributeError:
        print(
            '\nWARNING: No or limited FBC support limited! Please update your libSBML to the latest version.\n'
        )

    try:
        FBC3_VARIABLE_TYPES = {libsbml.FBC_VARIABLE_TYPE_LINEAR : 'linear',
                              libsbml.FBC_VARIABLE_TYPE_QUADRATIC : 'quadratic',
                              libsbml.FBC_VARIABLE_TYPE_INVALID : 'invalid',
                              'linear' : libsbml.FBC_VARIABLE_TYPE_LINEAR,
                              'quadratic' : libsbml.FBC_VARIABLE_TYPE_QUADRATIC,
                              'invalid' : libsbml.FBC_VARIABLE_TYPE_INVALID
                              }
    except AttributeError:
        print(
            '\nWARNING: Your libSBML does not contain FBCv3 support, please upgrade.'
        )




except ImportError:
    print(
        'WARNING: SBML support not available, please install libSBML, Python bindings with FBC (sbml.org)'
    )
    SBML_TYPES = {}
    _HAVE_SBML_ = False
try:
    GROUP_KINDS = {
        libsbml.GROUP_KIND_CLASSIFICATION: 'classification',
        libsbml.GROUP_KIND_PARTONOMY: 'partonomy',
        libsbml.GROUP_KIND_COLLECTION: 'collection',
        libsbml.GROUP_KIND_UNKNOWN: 'collection',
    }
    SBML_TYPES['group'] = libsbml.SBML_GROUPS_GROUP
except (AttributeError, NameError):
    print(
        'INFO: SBML+GROUPS support not available, update to latest version of libSBML if required'
    )
    GROUP_KINDS = {}
    _HAVE_GROUPS_ = False

_TEMP_XML_FILE_ = '_tmpxml.tmp'
FBA_NS = 'http://www.sbml.org/sbml/level3/version1/fba/version1'
METAPREFIX = 'meta_'

UNIT_DICTIONARY_L2 = {
    'area': {0: {'exponent': 2, 'kind': 'metre', 'multiplier': 1.0, 'scale': 0}},
    'length': {0: {'exponent': 1, 'kind': 'metre', 'multiplier': 1.0, 'scale': 0}},
    'substance': {0: {'exponent': 1, 'kind': 'mole', 'multiplier': 1.0, 'scale': 0}},
    'time': {0: {'exponent': 1, 'kind': 'second', 'multiplier': 1.0, 'scale': 0}},
    'volume': {0: {'exponent': 1, 'kind': 'litre', 'multiplier': 1.0, 'scale': 0}},
    'mmol_per_gDW_per_hr': {
        0: {'exponent': 1, 'kind': 'mole', 'multiplier': 1.0, 'scale': -3},
        1: {'exponent': -1, 'kind': 'gram', 'multiplier': 1.0, 'scale': 0},
        2: {'exponent': -1, 'kind': 'second', 'multiplier': 3600, 'scale': 0},
    },
}

UNIT_DICTIONARY = {
    'hour': {0: {'exponent': 1, 'kind': 'second', 'multiplier': 3600, 'scale': 0}},
    'mmol_per_gdw': {
        0: {'exponent': 1, 'kind': 'mole', 'multiplier': 1.0, 'scale': -3},
        1: {'exponent': -1, 'kind': 'gram', 'multiplier': 1.0, 'scale': 0},
    },
    'mmol_per_hour': {
        0: {'exponent': 1, 'kind': 'mole', 'multiplier': 1.0, 'scale': -3},
        1: {'exponent': -1, 'kind': 'second', 'multiplier': 3600, 'scale': 0},
    },
    'per_hour': {0: {'exponent': -1, 'kind': 'second', 'multiplier': 3600, 'scale': 0}},
}
MODEL_UNITS = {'extent': 'mmol_per_gdw', 'substance': 'mmol_per_gdw', 'time': 'hour'}

SBML_NS = [
    ('http://www.sbml.org/sbml/level3/version1/fbc/version3', 'L3V1FBC3'),
    ('http://www.sbml.org/sbml/level3/version1/fbc/version2', 'L3V1FBC2'),
    ('http://www.sbml.org/sbml/level3/version1/fbc/version1', 'L3V1FBC1'),
    ('http://www.sbml.org/sbml/level3/version2/fbc/version2', 'L3V2FBC2'),
    ('http://www.sbml.org/sbml/level3/version2/fbc/version1', 'L3V2FBC1'),
    ('http://www.sbml.org/sbml/level3/version2/core', 'L3V2core'),
    ('http://www.sbml.org/sbml/level3/version1/core', 'L3V1core'),
    ('http://www.sbml.org/sbml/level2/version4', 'L2'),
    ('http://www.sbml.org/sbml/level2', 'L2'),
]


"""
print libsbml.BQB_ENCODES            , 8  # "encodes",
print libsbml.BQB_HAS_PART           , 1  # "hasPart",
print libsbml.BQB_HAS_PROPERTY       , 10 # "hasProperty",
print libsbml.BQB_HAS_VERSION        , 4  # "hasVersion",
print libsbml.BQB_IS                 , 0  # "isA",
print libsbml.BQB_IS_DESCRIBED_BY    , 6  # "isDescribedBy",
print libsbml.BQB_IS_ENCODED_BY      , 7  # "isEncodedBy",
print libsbml.BQB_IS_HOMOLOG_TO      , 5  # "isHomologTo",
print libsbml.BQB_IS_PART_OF         , 2  # "isPartOf",
print libsbml.BQB_IS_PROPERTY_OF     , 11 # "isPropertyOf",
print libsbml.BQB_IS_VERSION_OF      , 3  # "isVersionOf",
print libsbml.BQB_OCCURS_IN          , 9  # "occursIn",
print libsbml.BQB_UNKNOWN            , 12 # None

print libsbml.BQM_IS                 , 0 # None
print libsbml.BQM_IS_DERIVED_FROM    , 2 # None
print libsbml.BQM_IS_DESCRIBED_BY    , 1 # None
print libsbml.BQM_UNKNOWN            , 3 # None
"""
BQ2CBMMAP = {
    8: "encodes",
    1: "hasPart",
    10: "hasProperty",
    4: "hasVersion",
    0: "is",
    6: "isDescribedBy",
    7: "isEncodedBy",
    5: "isHomologTo",
    2: "isPartOf",
    11: "isPropertyOf",
    3: "isVersionOf",
    9: "occursIn",
    12: None,
}

CBM2BQMAP = {
    "encodes": 8,
    "hasPart": 1,
    "hasProperty": 10,
    "hasVersion": 4,
    "is": 0,
    "isDescribedBy": 6,
    "isEncodedBy": 7,
    "isHomologTo": 5,
    "isPartOf": 2,
    "isPropertyOf": 11,
    "isVersionOf": 3,
    "occursIn": 9,
    None: 12,
}

BQM2CBMMAP = {0: "is", 1: "isDescribedBy", 2: "isDerivedFrom", 3: None}

CBM2BQMMAP = {"is": 0, "isDescribedBy": 1, "isDerivedFrom": 2, None: 3}

re_html_p = re.compile("<p>.*?</p>")
re_html_cobra_p = re.compile("<html:p>.*?</html:p>")
re_html_span = re.compile("<span>.*?</span>")


class MLStripper(HTMLParser):
    """
    Class for stripping a string of HTML/XML used from:
    http://stackoverflow.com/questions/753052/strip-html-from-strings-in-python

    """

    def __init__(self):
        """
        Strip a string of HTML/XML tags

        """
        HTMLParser.__init__(self)
        self.reset()
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        data = ''
        for d_ in self.fed:
            if d_.strip() != '':
                data += d_
        self.fed = []
        return data

# with Python 3.9+ moving unescape we now need to implement a nasty hack
__tagStripper__ = MLStripper()
try:
    __tagStripper__.unescape('dfsdfsdfsdfsdfsdfsdef')
except AttributeError as err:
    from html import unescape
    __tagStripper__.unescape = unescape


def xml_stripTags(html):
    """
    Strip a string of HTML/XML, returns a string

     - *html* the string containing html

    """
    __tagStripper__.feed(html)
    return __tagStripper__.get_data()


def formatSbmlId(s):
    out = ''
    for x in s:
        if x.isalnum():
            out += x
        else:
            out += '_'
    if not out[0].isalpha() and out[0] != '_':
        out = '_' + out
    return out


def sbml_readSBML2FBA(
    fname, work_dir=None, return_sbml_model=False, fake_boundary_species_search=False
):
    """
    Read in an SBML Level 2 file with FBA annotation where and return either a CBM model object
    or a (cbm_mod, sbml_mod) pair if return_sbml_model=True

     - *fname* is the filename
     - *work_dir* is the working directory (only used if not None)
     - *return_sbml_model* [default=False] return a a (cbm_mod, sbml_mod) pair
     - *fake_boundary_species_search* [default=False] after looking for the boundary_condition of a species search for overloaded id's <id>_b

    """

    assert (
        _HAVE_SBML_
    ), "\nSBML not available ... install libSBML with Python bindings for SBML support"
    if work_dir != None:
        fname = os.path.join(work_dir, fname)
    D = libsbml.readSBML(str(fname))
    M = D.getModel()
    if M == None:
        print('Invalid SBML, no model will be generated')
        if not return_sbml_model:
            return None
        else:
            return None, D
    F = csio.StringIO()
    model_id = M.getId()
    model_name = M.getName()
    if libsbml.getLibSBMLVersion() >= 51903:
        node_txt = libsbml.XMLNode.convertXMLNodeToString(M.getNotes())
    else:
        node_txt = libsbml.XMLNode_convertXMLNodeToString(M.getNotes())
    model_description = node_txt
    model_description = xml_stripTags(model_description).strip()

    # print(model_description)

    __HAVE_FBA_ANOT__ = False
    __HAVE_FBA_ANOT_OBJ__ = False
    __HAVE_FBA_ANOT_BNDS__ = False
    ANOT = M.getAnnotationString()
    if len(ANOT) > 0:

        __HAVE_FBA_ANOT__ = True
        if '<fba:listOfConstraints>' in ANOT:
            __HAVE_FBA_ANOT_BNDS__ = True
        # if '<fba:listOfObjectives/>' in ANOT:
        if '<fba:listOfFluxes>' in ANOT:
            __HAVE_FBA_ANOT_OBJ__ = True
        # fix some compatability stuff
        ANOT = ANOT.replace(
            '<annotation>',
            '<annotation xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/" xmlns:fba="http://www.sbml.org/sbml/level3/version1/fba/version1">',
        )
        # _TEMP_XML_FILE_ = createTempFileName()
        # _TEMP_XML_FILE_ = os.path.join(work_dir, _TEMP_XML_FILE_)
        # F = open(_TEMP_XML_FILE_,'w')
        F.write(ANOT)
        F.seek(0)
        # F.flush()
        # F.close()
    del ANOT
    if not __HAVE_FBA_ANOT_BNDS__ and not __HAVE_FBA_ANOT_OBJ__:
        __HAVE_FBA_ANOT__ = False

    SPEC = []
    for s in range(M.getNumSpecies()):
        SBSp = M.getSpecies(s)
        sid = SBSp.getId()
        boundCon = False
        if SBSp.getBoundaryCondition():
            if __DEBUG__:
                print('Real boundary metabolite: {}'.format(sid))
            boundCon = True
        CF = ''  # chemical formula
        CHRG = SBSp.getCharge()
        # if CHRG == 0:
        # CHRG = None
        NM = SBSp.getName()  # get name
        # to strip a BiGG file see CBTools
        S = CBModel.Species(
            SBSp.getId(),
            boundary=boundCon,
            name=NM,
            value=SBSp.getInitialConcentration(),
            compartment=SBSp.getCompartment(),
            charge=CHRG,
            chemFormula=CF,
        )
        # process notes field, get rid of <head>, <body> elements
        if libsbml.getLibSBMLVersion() >= 51903:
            specNotes = libsbml.XMLNode.convertXMLNodeToString(SBSp.getNotes())
        else:
            specNotes = libsbml.XMLNode_convertXMLNodeToString(SBSp.getNotes())
        S.annotation = sbml_readCOBRANote(specNotes)

        # Note: chemFormula works will have to see about charge GETFROMNAME!!!
        processSpeciesChargeChemFormulaAnnot(S, getFromName=True, overwriteCharge=True)

        # S.setAnnotation('note', specNotes)
        manot = sbml_getCVterms(SBSp, model=False)
        if manot != None:
            S.miriam = manot
        del manot
        SPEC.append(S)

    boundary_species = [s.getId() for s in SPEC if s.is_boundary]
    if len(boundary_species) == 0:
        if not fake_boundary_species_search:
            print(
                '\nINFO: No boundary species detected, if this is not what you expect try searching for boundary species using <name>_b (with \"fake_boundary_species_search=True\")?'
            )
        else:
            for s in range(M.getNumSpecies()):
                SBSp = M.getSpecies(s)
                sid = SBSp.getId()
                if sid[-2:] == '_b':
                    print('Fake boundary (_b) metabolite added: {}'.format(sid))
                    SPEC[s].is_boundary = True
    boundary_species = [s.getId() for s in SPEC if s.is_boundary]
    spec_id = [s.getId() for s in SPEC]

    REAC = []
    reactionIDs = []
    reactionsReversability = []
    for r in range(M.getNumReactions()):
        SBRe = M.getReaction(r)
        R_id = SBRe.getId()
        reagents = {}
        EXREAC = False
        reactionIDs.append(R_id)
        for rea in range(SBRe.getNumReactants()):
            spec = SBRe.getReactant(rea).getSpecies()
            stoi = -SBRe.getReactant(rea).getStoichiometry()
            # reagents.append((stoi,spec))
            if spec not in reagents:
                reagents[spec] = float(stoi)
            else:
                reagents[spec] += float(stoi)
            if spec in boundary_species:
                EXREAC = True
        for pr in range(SBRe.getNumProducts()):
            spec2 = SBRe.getProduct(pr).getSpecies()
            stoi2 = SBRe.getProduct(pr).getStoichiometry()
            # reagents.append((stoi2, spec2))
            if spec2 not in reagents:
                reagents[spec2] = float(stoi2)
            else:
                reagents[spec2] += float(stoi2)
            if spec2 in boundary_species:
                EXREAC = True
        R = CBModel.Reaction(
            SBRe.getId(), SBRe.getName(), reversible=SBRe.getReversible()
        )
        reactionsReversability.append(SBRe.getReversible())
        for r in reagents:
            R.addReagent(
                CBModel.Reagent('{}_{}'.format(SBRe.getId(), r), r, reagents[r])
            )
            # TODO check this; I'm almost sure it is nore needed anymore
            # if R.getId() not in SPEC[spec_id.index(r[1])].reagent_of:
            # SPEC[spec_id.index(r[1])].reagent_of.append(R.getId())
        del reagents
        if EXREAC:
            R.is_exchange = True

        if libsbml.getLibSBMLVersion() >= 51903:
            reacNotes = libsbml.XMLNode.convertXMLNodeToString(SBRe.getNotes())
        else:
            reacNotes = libsbml.XMLNode_convertXMLNodeToString(SBRe.getNotes())
        R.annotation = sbml_readCOBRANote(reacNotes)
        manot = sbml_getCVterms(SBRe, model=False)
        if manot != None:
            R.miriam = manot
        del manot
        REAC.append(R)

    CONSTR = []
    OBJFUNCout = []
    objfunc_data = {}
    multiobj = []
    if __HAVE_FBA_ANOT__:
        # root = ELTree.ElementTree(file=os.path.join(work_dir, _TEMP_XML_FILE_))
        root = ELTree.ElementTree(file=F)
        root_i = root.iter()
        for ri in root_i:
            if (
                ri.tag
                == '{http://www.sbml.org/sbml/level3/version1/fba/version1}fluxBalance'
            ):
                if __DEBUG__:
                    print(ri.tag)
                rootfba = ELTree.ElementTree(ri)
                root_fba_i = rootfba.iter()
        constraints = []
        for ret in root_fba_i:
            if (
                __HAVE_FBA_ANOT_BNDS__
                and ret.tag
                == '{http://www.sbml.org/sbml/level3/version1/fba/version1}listOfConstraints'
            ):
                if __DEBUG__:
                    print(ret.tag)
                for c in list(ret):
                    if __DEBUG__:
                        print('\t{}'.format(c.tag))
                    attrib = c.attrib
                    for a in attrib:
                        attrib.update(
                            {
                                a.replace(
                                    '{http://www.sbml.org/sbml/level3/version1/fba/version1}',
                                    '',
                                ): attrib.pop(a)
                            }
                        )
                    constraints.append(attrib)
            elif (
                __HAVE_FBA_ANOT_OBJ__
                and ret.tag
                == '{http://www.sbml.org/sbml/level3/version1/fba/version1}listOfObjectives'
            ):
                if __DEBUG__:
                    print(ret.tag)
                if len(ret) == 0:
                    print('No objectives in listOfObjectives')
                else:
                    activeId = None
                    if (
                        '{http://www.sbml.org/sbml/level3/version1/fba/version1}activeObjective'
                        in ret.attrib
                    ):
                        activeId = ret.attrib[
                            '{http://www.sbml.org/sbml/level3/version1/fba/version1}activeObjective'
                        ]
                    for obj in list(ret):
                        ##  print obj.attrib
                        if (
                            '{http://www.sbml.org/sbml/level3/version1/fba/version1}type'
                            in obj.attrib
                        ):
                            ftype = obj.attrib[
                                '{http://www.sbml.org/sbml/level3/version1/fba/version1}type'
                            ]
                        else:
                            ftype = None
                        ##  raw_input(ftype)
                        sid = obj.attrib['id']
                        # multiobj = []
                        for c_ in list(obj):
                            fo = []
                            for cc_ in list(c_):
                                fo.append(cc_.attrib)
                            multiobj.append(fo)

                        if __DEBUG__:
                            print(objfunc_data)
                        obj_ = multiobj[-1]
                        for flobj_ in obj_:
                            for a in tuple(flobj_):
                                flobj_.update(
                                    {
                                        a.replace(
                                            '{http://www.sbml.org/sbml/level3/version1/fba/version1}',
                                            '',
                                        ): flobj_.pop(a)
                                    }
                                )
                            flobj_.update({'type': ftype})
                            flobj_.update({'id': sid})

        cntr = 0
        boundReactionIDs = []
        LboundReactionIDs = []
        UboundReactionIDs = []
        AboundReactionIDs = []
        DefinedReactionIDs = []
        for c in constraints:
            if 'id' in c:
                newId = c['id']
                DefinedReactionIDs.append(c['reaction'])
            else:
                O = 'C_%i' % cntr
                if c['operation'] in ['greater', 'greaterEqual', '>', '>=']:
                    O = 'lower'
                    LboundReactionIDs.append(c['reaction'])
                elif c['operation'] in ['less', 'lessEqual', '<', '<=']:
                    O = 'upper'
                    UboundReactionIDs.append(c['reaction'])
                elif c['operation'] in ['equal', '=']:
                    O = 'equal'
                    AboundReactionIDs.append(c['reaction'])
                newId = '%s_%s_bnd' % (c['reaction'], O)
            CONSTR.append(
                CBModel.FluxBound(
                    newId, c['reaction'], c['operation'], float(c['value'])
                )
            )
            cntr += 1
            if c['reaction'] not in boundReactionIDs:
                boundReactionIDs.append(c['reaction'])

        # undefined flux bounds are given infinite value
        ubcntr = 0
        ##  print boundReactionIDs
        ##  print reactionIDs

        for J in range(len(reactionIDs)):
            ##  print reactionIDs[J], reactionsReversability[J]
            LBt = False
            UBt = False
            ABt = False
            if reactionIDs[J] in DefinedReactionIDs:
                LBt = UBt = ABt = True
            if reactionIDs[J] in LboundReactionIDs:
                LBt = True
            if reactionIDs[J] in UboundReactionIDs:
                UBt = True
            if reactionIDs[J] in AboundReactionIDs:
                ABt = True
            if not LBt and not UBt and not ABt:
                # print 'not LBt and not UBt and not ABt'
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'lower')
                if reactionsReversability[J]:
                    ##  print 'Adding reversible'
                    CONSTR.append(
                        CBModel.FluxBound(
                            newId, reactionIDs[J], 'greaterEqual', -numpy.inf
                        )
                    )
                else:
                    ##  print 'Adding irreversible'
                    CONSTR.append(
                        CBModel.FluxBound(newId, reactionIDs[J], 'greaterEqual', 0.0)
                    )
                ubcntr += 1
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'upper')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'lessEqual', numpy.inf)
                )
                ubcntr += 1
            elif not LBt and not ABt:
                # print 'not LBt and not ABt'
                # print reactionIDs[J]
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'lower')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'greaterEqual', -numpy.inf)
                )
                ubcntr += 1
                # print 'Added new lower bound', newId
            elif not UBt and not ABt:
                # print 'not UBt and not ABt'
                # print reactionIDs[J]
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'upper')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'lessEqual', numpy.inf)
                )
                ubcntr += 1
                # print 'Added new upper bound', newId

    # printl(CONSTR)
    # printl(__HAVE_FBA_ANOT_BNDS__)
    if not __HAVE_FBA_ANOT_BNDS__ and len(CONSTR) == 0:
        for r_ in REAC:
            rid = r_.getId()
            if r_.reversible:
                newId = '{}_{}_bnd'.format(rid, 'lower')
                CONSTR.append(CBModel.FluxBound(newId, rid, 'greaterEqual', -numpy.inf))
                newId = '{}_{}_bnd'.format(r_.getId(), 'upper')
                CONSTR.append(CBModel.FluxBound(newId, rid, 'lessEqual', numpy.inf))
            else:
                newId = '{}_{}_bnd'.format(rid, 'lower')
                CONSTR.append(CBModel.FluxBound(newId, rid, 'greaterEqual', 0.0))
                newId = '{}_{}_bnd'.format(r_.getId(), 'upper')
                CONSTR.append(CBModel.FluxBound(newId, rid, 'lessEqual', numpy.inf))

    # build model
    fm = CBModel.Model(model_id)
    fm.sourcefile = fname
    manot = sbml_getCVterms(M, model=True)
    if manot != None:
        fm.miriam = manot
    del manot
    sbmh = M.getModelHistory()
    if sbmh != None:
        cdate = sbmh.getCreatedDate()
        if cdate != None:
            cdate = (
                cdate.getYear(),
                cdate.getMonth(),
                cdate.getDay(),
                cdate.getHour(),
                cdate.getMinute(),
                cdate.getSecond(),
            )
            fm.setCreatedDate(cdate)
        mdate = sbmh.getModifiedDate()
        if mdate != None:
            mdate = (
                mdate.getYear(),
                mdate.getMonth(),
                mdate.getDay(),
                mdate.getHour(),
                mdate.getMinute(),
                mdate.getSecond(),
            )
            fm.setModifiedDate(mdate)
        for m_ in range(sbmh.getNumCreators()):
            sbc = sbmh.getCreator(m_)
            fm.addModelCreator(
                sbc.getGivenName(),
                sbc.getFamilyName(),
                sbc.getOrganisation(),
                sbc.getEmail(),
            )
    fm.name = model_name
    fm.description = model_description
    for s in SPEC:
        fm.addSpecies(s)
    for r in REAC:
        fm.addReaction(r, create_default_bounds=False)
    for c in CONSTR:
        fm.addFluxBound(c)
    if len(multiobj) > 0:
        otype = None
        oid = None
        for obj_ in multiobj:
            flobjs = []
            for fobj_ in obj_:
                otype = fobj_['type']
                oid = fobj_['id']
                flobjs.append(
                    CBModel.FluxObjective(
                        '{}_{}_flobj'.format(fobj_['id'], fobj_['reaction']),
                        fobj_['reaction'],
                        float(fobj_['coefficient']),
                    )
                )
            o = CBModel.Objective(oid, otype)
            if len(multiobj) == 1:
                fm.addObjective(o, active=True)
            elif o.getId() == activeId:
                fm.addObjective(o, active=True)
            else:
                fm.addObjective(o, active=False)
            OBJFUNCout.append(o)
            print(o.getId(), activeId)
            for f_ in flobjs:
                o.addFluxObjective(f_)
    del SPEC, REAC, CONSTR
    try:
        F.close()
        del F
    except Exception as ex:
        print(ex)
        print('INFO: File deletion failure')
    if not __HAVE_FBA_ANOT__:
        print('\nWARNING: Missing FBA annotations found')
        if not __HAVE_FBA_ANOT_BNDS__:
            print('WARNING: No FluxBounds defined (generic bounds created)')
        if not __HAVE_FBA_ANOT_OBJ__:
            print('WARNING: No Objective function defined please create one')
        time.sleep(1)
    fm._SBML_LEVEL_ = 2
    try:
        fm.buildStoichMatrix()
    except Exception as ex:
        print(ex)
        print('INFO: unable to construct stoichiometric matrix')
    if not return_sbml_model:
        del M, D
        return fm
    else:
        return fm, M.clone()


def xml_createSBML2FBADoc():
    """
    Create a 'document' to store the SBML2FBA annotation, returns:

     - *DOC* a minidom document

    """
    DOM = getDOMImplementation()
    DOC = DOM.createDocument(FBA_NS, 'fba:fluxBalance', None)
    DOC.documentElement.setAttribute('xmlns:fba', FBA_NS)
    TOP = DOC.documentElement
    LoC = DOC.createElementNS(FBA_NS, 'fba:listOfConstraints')
    TOP.appendChild(LoC)
    LoO = DOC.createElementNS(FBA_NS, 'fba:listOfObjectives')
    TOP.appendChild(LoO)
    return DOC


def xml_viewSBML2FBAXML(document, fname=None):
    """
    Print a minidom XML document to screen or file, arguments:

     - *document* a minidom XML document
     - *fname* [default=None] by default print to screen or write to file fname

    """
    ##  print DOC.toxml('UTF-8')
    ##  print document.toprettyxml(indent=' ',newl='\n',encoding='UTF-8')
    if fname != None:
        FO = open(fname, 'w')
        document.writexml(FO, indent='', addindent='    ', newl='\n', encoding='UTF-8')
        FO.flush()
        FO.close()


def xml_addSBML2FBAFluxBound(document, rid, operator, value, fbid=None):
    """
    Adds an SBML3FBA flux bound to the document:

     - *document* a minidom XML document created by xml_createSBML2FBADoc
     - *rid* the reaction id
     - *operator* one of ['greater','greaterEqual','less','lessEqual','equal','>','>=','<','<=','=']
     - *value* a float which will be cast to a string using str(value)
     - *fbid* the flux bound id, autogenerated by default

    """
    OPER = [
        'greater',
        'greaterEqual',
        'less',
        'lessEqual',
        'equal',
        '>',
        '>=',
        '<',
        '<=',
        '=',
    ]
    assert operator in OPER, '\nBad operator %s should be one of %s' % (
        operator,
        str(OPER),
    )
    LoC = document.getElementsByTagName('fba:listOfConstraints')[0]
    F = document.createElementNS(FBA_NS, 'fba:constraint')
    F.setAttributeNS(FBA_NS, 'fba:reaction', rid)
    F.setAttributeNS(FBA_NS, 'fba:operation', operator)
    value = float(value)
    if numpy.isinf(value):
        if value > 0:
            value = 'inf'
        else:
            value = '-inf'
    else:
        value = str(value)
    F.setAttributeNS(FBA_NS, 'fba:value', value)
    if fbid == None:
        O = 'type'
        if operator in ['greater', 'greaterEqual', '>', '>=']:
            O = 'lower'
        elif operator in ['less', 'lessEqual', '<', '<=']:
            O = 'upper'
        elif operator in ['equal', '=']:
            O = 'equal'
        fbid = '%s_%s_bnd' % (rid, O)

    # METAID
    # F.setAttributeNS(FBA_NS,'metaid', fbid)
    F.setAttributeNS(FBA_NS, 'id', fbid)
    LoC.appendChild(F)


def xml_createListOfFluxObjectives(document, l):
    """
    Create a list of fluxObjectives to add to an Objective:

     - *document* a minidom XML document created by xml_createSBML2FBADoc
     - *fluxobjs* a list of (rid, coefficient) tuples

    """
    LoF = document.createElementNS(FBA_NS, 'fba:listOfFluxes')
    for F in fluxObjectives:
        FO = document.createElementNS(FBA_NS, 'fba:fluxObjective')
        FO.setAttributeNS(FBA_NS, 'fba:reaction', F[0])
        FO.setAttributeNS(FBA_NS, 'fba:coefficient', F[1])
        LoF.appendChild(FO)
    return LoF


def xml_createSBML2FBAObjective(document, oid, sense, fluxObjectives):
    """
    Create a list of fluxObjectives to add to an Objective:

     - *document* a minidom XML document created by xml_createSBML2FBADoc
     - *oid* the objective id
     - *sense* a string containing the objective sense either: **maximize** or **minimize**
     - *fluxObjectives* a list of (rid, coefficient) tuples

    """
    sense = sense.lower()
    if sense == 'max':
        sense = 'maximize'
    if sense == 'min':
        sense = 'minimize'
    if sense in ['maximise', 'minimise']:
        sense = sense.replace('se', 'ze')
    assert sense in ['maximize', 'minimize'], (
        "\nType must be ['maximize', 'minimize'] not %s" % sense
    )
    OBJ = document.createElementNS(FBA_NS, 'fba:objective')
    OBJ.setAttributeNS(FBA_NS, 'id', oid)
    # METAID
    # OBJ.setAttributeNS(FBA_NS,'metaid', oid)
    OBJ.setAttributeNS(FBA_NS, 'fba:type', sense)
    LoF = document.createElementNS(FBA_NS, 'fba:listOfFluxes')
    for F in fluxObjectives:
        FO = document.createElementNS(FBA_NS, 'fba:fluxObjective')
        FO.setAttributeNS(FBA_NS, 'fba:reaction', F[0])
        FO.setAttributeNS(FBA_NS, 'fba:coefficient', str(F[1]))
        LoF.appendChild(FO)
    OBJ.appendChild(LoF)
    return OBJ


def xml_addSBML2FBAObjective(document, objective, active=True):
    """
    Adds an objective element to the documents listOfObjectives and sets the active attribute:

     - *document* a minidom XML document created by `xml_createSBML2FBADoc`
     - *objective* a minidom XML objective element created with `xml_createSBML2FBAObjective`
     - *active* [default=True] a boolean flag specifiying whether this objective is active

    """

    LoO = document.getElementsByTagName('fba:listOfObjectives')[0]
    if active:
        xid = objective.getAttributeNS(FBA_NS, 'id')
        LoO.setAttributeNS(FBA_NS, 'fba:activeObjective', xid)
        print('Setting active objective: {}'.format(xid))
    LoO.appendChild(objective)


def xml_getSBML2FBAannotation(fba, fname=None):
    """
    Takes an FBA model object and returns the SBML3FBA annotation as an XML string:

     - *fba* an fba model object
     - *fname* [default=None] if supplied the XML will be written to file *fname*

    """

    DOC = xml_createSBML2FBADoc()
    for f in fba.flux_bounds:
        ##  print f.value
        xml_addSBML2FBAFluxBound(DOC, f.reaction, f.operation, f.value, f.getId())

    for o in fba.objectives:
        fluxobjs = [(fo.reaction, fo.coefficient) for fo in o.flux_objectives]
        OBJ = xml_createSBML2FBAObjective(DOC, o.id, o.operation, fluxobjs)
        if o.id == fba.objectives[fba.activeObjIdx].id:
            xml_addSBML2FBAObjective(DOC, OBJ, active=True)
        else:
            xml_addSBML2FBAObjective(DOC, OBJ, active=False)
    if fname != None:
        xml_viewSBML2FBAXML(DOC, fname)

    return DOC.toprettyxml(indent=' ', newl='\n', encoding='UTF-8')


def sbml_createModelL2(fba, level=2, version=1):
    """
    Create an SBML model and document:

     - *fba* a PySCeSCBM model instance
     - *level* always 2
     - *version* always 1

    and returns:

     - *model* an SBML model

    """
    SBML_LEVEL = level = 2
    SBML_VERSION = version = 1

    if fba.getId() == '' or fba.getId() == None:
        mid0 = 'FBAModel'
    else:
        mid0 = fba.getId()

    mid = ''
    for l in mid0:
        if l.isalnum():
            mid += l
        else:
            mid += '_'

    if not mid[0].isalpha():
        mid = 'id_' + mid

    mname = fba.getName()
    if mname in [None, '', ' ']:
        mname = 'SBML_CB_MODEL'

    if libsbml.getLibSBMLVersion() < 40000:
        document = libsbml.SBMLDocument()
        document.setLevelAndVersion(SBML_LEVEL, SBML_VERSION)
    else:
        # document = libsbml.SBMLDocument(SBML_LEVEL, SBML_VERSION)
        document = libsbml.SBMLDocument(SBML_LEVEL, SBML_VERSION)
    document.getNamespaces().add("http://www.w3.org/1999/xhtml", "html")

    model = document.createModel(fba.getId())
    model.setMetaId(METAPREFIX + fba.getId())
    ## can't do this right now with the custom annotations
    # miriam = fba.miriam.getAllMIRIAMUris()
    # print 'miriam', miriam
    # if len(miriam) > 0:
    # sbml_setCVterms(model, miriam, model=True)
    ## while this may look right it totally screws up things later, notably reaction reagents don't get set
    # ns = model.getNamespaces()
    # ns.add("http://www.w3.org/1999/xhtml", "html")
    # model.setNamespaces(ns)
    # print model.getNamespaces()

    ##  document = libsbml.SBMLDocument(SBML_LEVEL, SBML_VERSION)
    return model, document


def sbml_setCompartmentsL2(model, compartments):
    if len(compartments) < 1:
        print('Warning: no compartments defined adding one called \"cell\"')
        compartments.update({'cell': {'id': 'cell', 'dimensions': 3, 'size': 1.0}})
    for cs in compartments:
        comp_def = model.createCompartment()
        comp_def.setId(compartments[cs]['id'])
        comp_def.setName(compartments[cs]['id'])
        comp_def.setVolume(float(compartments[cs]['size']))
        comp_def.setConstant(False)


def sbml_setCompartmentsL3(model, fba):
    '''
    Sets the model compartments.

     - *model* a libSBML model instance
     - *fba* a PySCeSCBM model instance

    '''
    for cs in fba.compartments:
        comp_def = model.createCompartment()
        comp_def.setId(cs.getId())
        comp_def.setName(cs.getName())
        # TODO: look into this for user selectable units
        comp_def.setUnits('dimensionless')

        size = cs.getSize()
        if size == None or size == '' or numpy.isnan(size):
            # TODO: need to decide what to do here
            comp_def.setSize(1)
        else:
            comp_def.setSize(cs.getSize())

        comp_def.setSpatialDimensions(cs.getDimensions())
        comp_def.setConstant(True)
        comp_def.setMetaId(METAPREFIX + cs.getId())

        if len(cs.getAnnotations()) > 0:
            annoSTRnew = sbml_writeKeyValueDataAnnotation(cs.getAnnotations())
            annores = comp_def.appendAnnotation(annoSTRnew)
            if annores == -3:
                print('Invalid annotation in reaction:', cs.getId())
                print(cs.getAnnotations())
        if cs.miriam != None:
            miriam = cs.miriam.getAllMIRIAMUris()
            if len(miriam) > 0:
                sbml_setCVterms(comp_def, miriam, model=False)
        sboterm = cs.getSBOterm()
        if sboterm is not None and sboterm != '':
            comp_def.setSBOTerm(str(sboterm))
        notes = cs.getNotes()
        if notes != '' and notes is not None:
            sbml_setNotes3(comp_def, notes)


def sbml_setParametersL3Fbc(fbcmod, add_cbmpy_anno=True, fbc_version=2):
    """
    Add non fluxbound related parameters to the model

     - *fbcmod* a CBM2SBML instance
     - *add_cbmpy_anno* [default=True] add CBMPy KeyValueData annotation.

    """

    cntr = 0
    for par in fbcmod.fba.parameters:
        if not par._is_fluxbound_:
            fbcmod.createParParameter(par, add_cbmpy_anno)
            cntr += 1

    print('INFO: added {} non fluxbound parameters to model'.format(cntr))


def sbml_setAnnotationsL3Fbc(cbmo, sbmlo):
    """
    Add CBMPy Fbase annotations to an SBML object, MIRIAM, SBO, Notes. Should
    be called last when creating SBML objects.

     - *cbmo* the CBMPy object
     - *sbmlo* SBML object

     Note: this function should be used for new code, old code still needs to be
     refactored.

    """

    if len(cbmo.getAnnotations()) > 0:
        annoSTRnew = sbml_writeKeyValueDataAnnotation(cbmo.getAnnotations())
        annores = sbmlo.appendAnnotation(annoSTRnew)
    if cbmo.miriam != None:
        miriam = cbmo.miriam.getAllMIRIAMUris()
        if len(miriam) > 0:
            sbml_setCVterms(sbmlo, miriam, model=False)
    sboterm = cbmo.getSBOterm()
    if sboterm is not None and sboterm != '':
        sbmlo.setSBOTerm(str(sboterm))
    notes = cbmo.getNotes()
    if notes != '' and notes is not None:
        sbml_setNotes3(sbmlo, notes)


def sbml_setDescription(model, fba):
    '''
    Sets the model description as a <note> containing `txt` in an HTML paragraph on the model object.

     - *model* a libSBML model instance
     - *fba* a PySCeSCBM model instance

    '''
    ##  try: UseR = getuser()
    ##  except: UseR = ''
    notes = ''
    if fba.getNotes().strip() in ['', None, ' ']:
        # notes = notes.replace('<notes>', '').replace('</notes>','')
        notes = '<p xmlns="http://www.w3.org/1999/xhtml"><br/><span size="small">Model \"<strong>{}</strong>\" ({}) generated with <a href="http://cbmpy.sourceforge.net">CBMPy</a> ({}) on {}.</span></p>'.format(
            fba.getId(),
            fba.getName(),
            __version__,
            time.strftime("%a, %d %b %Y %H:%M:%S"),
        )
    else:
        notes = fba.getNotes().strip()
        # notes = '<p xmlns="http://www.w3.org/1999/xhtml"><span style="font-family: Courier New,Courier,monospace;">{}</span></p>\n'.format(fba.getNotes().strip())
    sbml_setNotes3(model, notes)

    if len(fba.getAnnotations()) > 0:
        annoSTRnew = sbml_writeKeyValueDataAnnotation(fba.getAnnotations())
        annores = model.appendAnnotation(annoSTRnew)
        if annores == -3:
            print('Invalid annotation in model:', fba.getId())
            print(fba.getAnnotations(), '\n')


def sbml_setNotes3(obj, s):
    """
    Formats the CBMPy notes as an SBML note and adds it to the SBMl object

     - *obj* an SBML object
     - *s* a string that should be added as a note

    """
    s = s.replace('<notes>', '').replace('</notes>', '')
    if '</html:body>' not in s and '</body>' not in s:
        s = '<body xmlns="http://www.w3.org/1999/xhtml">{}</body>'.format(s)

    # this makes sure that whatever is there is encoded into ASCII
    # this is a workaround to #15 until I work out wtf is going on
    if os.sys.version_info[0] == 3:  # and libsbml.LIBSBML_VERSION == 51700:
        s = s.encode(encoding='ascii', errors='replace')
        s = s.decode('utf-8', errors='ignore')

    res = obj.setNotes(s)
    if res != 0:
        print(res, s)
        return False
    else:
        return True


def sbml_getNotes(obj):
    """
    Returns the SBML objects notes

    - *obj* an SBML object

    """
    notes = ''
    try:
        if libsbml.getLibSBMLVersion() >= 51903:
            notes = libsbml.XMLNode.convertXMLNodeToString(obj.getNotes())
        else:
            notes = libsbml.XMLNode_convertXMLNodeToString(obj.getNotes())
        if notes != '' and notes is not None:
            # too aggressive but efficient behaviour removed
            # notes = xml_stripTags(notes).strip()
            # notes = notes.replace('<html:body>', '').replace('</html:body>','')cmod.getDescription().replace('\xa0', ' ')
            notes = notes.replace('<notes>', '').replace('</notes>', '')
            # this removes urlencoded &nbsp; in xhtml (\xa0 don't ^%&*() ask!!!) which mistranslates when mapped through python unicode
            # notes = notes.replace('\xa0', ' ')
    except Exception as why:
        print(why)
    return notes


def sbml_setUnits(model, units=None, give_default=False, L3=True):
    """
    Adds units to the model:

    - *model* a libSBML model instance
    - *units* [default=None] a dictionary of units, if None default units are used
    - *give_default* [default=False] if true method returns the default unit dictionary
    - *L3* [default=True] use the L3 defaults

    """

    if units is None:
        if L3:
            units = UNIT_DICTIONARY
        else:
            units = UNIT_DICTIONARY_L2

    if give_default:
        return units

    for un in units:
        vdef = model.createUnitDefinition()
        vdef.setId(un)
        vdef.setName(un)
        if L3 and un in MODEL_UNITS.values():
            if MODEL_UNITS['time'] == un:
                model.setTimeUnits(un)
            if MODEL_UNITS['extent'] == un:
                model.setExtentUnits(un)
            if MODEL_UNITS['substance'] == un:
                model.setSubstanceUnits(un)
        for u in range(len(units[un])):
            vu = model.createUnit()
            vu.setKind(libsbml.UnitKind_forName(units[un][u]['kind']))
            vu.setMultiplier(units[un][u]['multiplier'])
            vu.setScale(int(units[un][u]['scale']))
            vu.setExponent(int(units[un][u]['exponent']))
            vu.setOffset(0)


def sbml_writeAnnotationsAsCOBRANote(annotations):
    """
    Writes the annotations dictionary as a COBRA compatible SBML <note>

    """
    annoSTR = ''
    for K in annotations:
        # this is to keep COBRA happy
        if K == 'subsystem':
            Kcob = 'SUBSYSTEM'
        elif K in ['gene_association', 'GENE_ASSOCIATION']:
            Kcob = 'GENE ASSOCIATION'
        else:
            Kcob = K
        val = str(annotations[K]).replace('<', '&lt;').replace('>', '&gt;')
        annoSTR += "<html:p>{}: {}</html:p>\n".format(Kcob, val)
    annoSTR = '{}'.format(annoSTR)
    return annoSTR


def sbml_setSpeciesL2(model, fba, return_dicts=False):
    """
    Add the species definitions to the SBML object:

     - *model* [default=''] a libSBML model instance or can be None if *return_dicts* == True
     - *fba* a PySCeSCBM model instance
     - *return_dicts* [default=False] only returns the compartment and species dictionaries without updated the SBML

    returns:

     - *compartments* a dictionary of compartments (except when give *return_dicts* argument)

    """

    compartments = {}
    species = {}
    for s in fba.species:
        if s.compartment not in compartments:
            if s.compartment == None or s.compartment == '':
                pass
            else:
                compartments.update(
                    {s.compartment: {'id': s.compartment, 'size': 1.0, 'dimensions': 3}}
                )
        # print '++', s.charge
        if s.miriam == None:
            miriam = {}
        else:
            miriam = s.miriam.getAllMIRIAMUris()
        species.update(
            {
                s.getId(): {
                    'id': s.getId(),
                    'compartment': s.compartment,
                    'name': s.getName(),
                    'charge': s.charge,
                    'value': s.value,
                    'annotation': s.getAnnotations().copy(),
                    'boundary': s.is_boundary,
                    'chemFormula': s.chemFormula,
                    'miriam': miriam,
                }
            }
        )

    if len(compartments) < 1:
        print('Warning: no compartments defined adding one called \"cell\"')
        compartments.update({'cell': {'id': 'cell', 'dimensions': 3, 'size': 1.0}})

    if return_dicts:
        return compartments, species

    keys = sorted(species)
    # keys.sort()
    for spe in keys:
        s = model.createSpecies()
        s.setId(species[spe]['id'])
        # METAID
        s.setMetaId(METAPREFIX + species[spe]['id'])
        s.setName(species[spe]['name'])
        # if not (species[spe]['charge'] != None or species[spe]['charge'] != ''):
        # s.setCharge(int(species[spe]['charge']))
        Ranno = species[spe]['annotation']
        if species[spe]['charge'] not in ['', None]:
            # print species[spe]['charge'], int(species[spe]['charge'])
            if (
                s.setCharge(int(species[spe]['charge']))
                != libsbml.LIBSBML_OPERATION_SUCCESS
            ):
                print(
                    'Unable to set charge, L{}V{} > L2V1'.format(
                        s.getLevel(), s.getVersion()
                    )
                )
                if species[spe]['charge'] not in ['', None, 0]:
                    Ranno.update({'CHARGE': species[spe]['charge']})

        if species[spe]['compartment'] == None or species[spe]['compartment'] == '':
            s.setCompartment('cell')
        else:
            s.setCompartment(species[spe]['compartment'])
        if species[spe]['boundary']:
            s.setBoundaryCondition(True)
            # s.setConstant(True)
        else:
            s.setBoundaryCondition(False)
        if species[spe]['value'] not in ['', None]:
            s.setInitialConcentration(float(species[spe]['value']))
            ##  s.setInitialAmount(float(species[spe]['value']))
        s.setHasOnlySubstanceUnits(True)

        annoSTR = None
        # if species[spe]['chemFormula'] != None and species[spe]['chemFormula'] != '':
        if species[spe]['chemFormula'] not in ['', None]:
            Ranno.update({'FORMULA': species[spe]['chemFormula']})
        if len(Ranno) > 0:
            annoSTR = sbml_writeAnnotationsAsCOBRANote(Ranno)
        if annoSTR != None:
            nres = s.setNotes(annoSTR)
        # this should go last
        if len(species[spe]['miriam']) > 0:
            sbml_setCVterms(s, species[spe]['miriam'], model=False)
    return compartments


def sbml_setReactionsL2(model, fba, return_dict=False):
    """
    Add the FBA instance reactions to the SBML model

     - *model* an SBML model instance
     - *fba* a PySCeSCBM model instance
     - *return_dict* [default=False] if True do not add reactions to SBML document instead return a dictionary description of the reactions

    """

    SBML_LEVEL = 2
    SBML_VERSION = 1
    reactions = {}

    for r in fba.reactions:
        reactants = []
        products = []
        for rg in r.getStoichiometry():
            if rg[0] > 0.0:
                products.append(rg)
            else:
                reactants.append(rg)
        if r.miriam == None:
            miriam = {}
        else:
            miriam = r.miriam.getAllMIRIAMUris()
        reactions.update(
            {
                r.getId(): {
                    'id': r.getId(),
                    'reactants': reactants,
                    'products': products,
                    'name': r.getName(),
                    'reversible': r.reversible,
                    'exchange': r.is_exchange,
                    'annotation': r.getAnnotations().copy(),
                    'compartment': r.compartment,
                    'miriam': miriam,
                }
            }
        )
    if return_dict:
        return reactions

    for rxn in reactions:
        # print 'Adding reaction:', reactions[rxn]['id']
        r = model.createReaction()
        r.setId(reactions[rxn]['id'])
        # METAID
        r.setMetaId(METAPREFIX + reactions[rxn]['id'])
        r.setName(reactions[rxn]['name'])
        for s in range(len(reactions[rxn]['reactants'])):
            # print '\t' + reactions[rxn]['id'] +' has substrate: ' + s.name + ' (%s)' % abs(rxn.stoichiometry[s.name])
            if libsbml.getLibSBMLVersion() < 40000:
                sref = libsbml.SpeciesReference(
                    reactions[rxn]['reactants'][s][1],
                    abs(float(reactions[rxn]['reactants'][s][0])),
                )
            else:
                sref = libsbml.SpeciesReference(SBML_LEVEL, SBML_VERSION)
                sref.setStoichiometry(abs(float(reactions[rxn]['reactants'][s][0])))
                sref.setSpecies(reactions[rxn]['reactants'][s][1])
            # METAID
            # sref.setMetaId('%s_%s' % (reactions[rxn]['id'], reactions[rxn]['reactants'][s][1]))
            assert r.addReactant(sref) == 0, '\nError setting reagent'
        for p in range(len(reactions[rxn]['products'])):
            # print '\t' + reactions[rxn]['id'] +' has product: ' + p.name + ' (%s)' % abs(rxn.stoichiometry[p.name])
            if libsbml.getLibSBMLVersion() < 40000:
                pref = libsbml.SpeciesReference(
                    reactions[rxn]['products'][p][1],
                    abs(float(reactions[rxn]['products'][p][0])),
                )
            else:
                pref = libsbml.SpeciesReference(SBML_LEVEL, SBML_VERSION)
                pref.setStoichiometry(abs(float(reactions[rxn]['products'][p][0])))
                pref.setSpecies(reactions[rxn]['products'][p][1])
            # METAID
            # pref.setMetaId('%s_%s' % (reactions[rxn]['id'], reactions[rxn]['products'][p][1]))
            assert r.addProduct(pref) == 0, '\nError setting product'
        Ranno = reactions[rxn]['annotation']
        # TODO: investigate this
        # if 'GENE_ASSOCIATION' in Ranno:
        # Ranno.update({'GENE ASSOCIATION' : Ranno.pop('GENE_ASSOCIATION')})
        # print 'Ranno', Ranno
        annoSTR = None
        if len(Ranno) > 0:
            annoSTR = sbml_writeAnnotationsAsCOBRANote(Ranno)
        if annoSTR != None:
            nres = r.setNotes(annoSTR)
            if nres != 0:
                print(nres, annoSTR)
            annoSTRnew = sbml_writeKeyValueDataAnnotation(
                Ranno
            )  # future compatability with L3
            annores = r.appendAnnotation(annoSTRnew)
        if reactions[rxn]['reversible']:
            r.setReversible(True)
        else:
            r.setReversible(False)
        # this should go last
        # this may be a bug in libSBML in that setCVterm needs to come after appendAnnotation
        if len(reactions[rxn]['miriam']) > 0:
            sbml_setCVterms(r, reactions[rxn]['miriam'], model=False)


def sbml_exportSBML2FBAModel(
    document, filename, directory=None, return_doc=False, remove_note_body=False
):
    """
    Writes an SBML model object to file. Note this is an internal SBML method use `sbml_writeSBML2FBA()` to write an FBA model:

     - *model* a libSBML model instance
     - *filename* the output filename
     - *directory* [default=None] by default use filename otherwise join, <dir><filename>
     - *return_doc* [default=False] return the SBML document used to write the XML

    """
    if directory != None:
        assert os.path.exists(directory), '\n%s does not exist.' % directory
        filename = os.path.join(directory, filename)
    print('Writing file: {}'.format(filename))

    try:
        UseR = getuser()
    except:
        UseR = ''

    SBML_LEVEL = 2
    SBML_VERSION = 1

    # raw_input('L%sV%s' % (document.getLevel(),document.getVersion()))
    h1 = '<?xml version="1.0" encoding="utf-8"?>\n'
    h1 += (
        '<!-- Created with CBMPy ('
        + __version__
        + ') on '
        + time.strftime("%a, %d %b %Y %H:%M:%S")
        + ' by '
        + UseR
        + ' -->\n'
    )

    SBMLDoc = document.toSBML()
    if remove_note_body:
        SBMLDoc = SBMLDoc.replace(
            '<body xmlns:html="http://www.w3.org/1999/xhtml">', ''
        )
        SBMLDoc = SBMLDoc.replace('</body>', '')

    F = open(filename, 'w')
    F.write(h1 + SBMLDoc)
    F.flush()
    F.close()
    print('Model exported as: {}'.format(filename))


def sbml_writeSBML2FBA(fba, fname, directory=None, sbml_level_version=None):
    """
    Takes an FBA model object and writes it to file as SBML L3 FBA:

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname
     - *sbml_level_version* [default=None] a tuple containing the SBML level and version e.g. (2,4) (ignored)

    """
    ##  __DEBUG__ = True

    assert (
        _HAVE_SBML_
    ), "\nSBML not available ... install libSBML with Python bindings for SBML support"

    SBML_LEVEL = 2
    SBML_VERSION = 1
    #  create L3 as annotation
    L3FBA = xml_getSBML2FBAannotation(fba, 'test_xml.xml')  # should be None

    # create a model
    SMOD, document = sbml_createModelL2(fba, level=SBML_LEVEL, version=SBML_VERSION)
    sbml_setDescription(SMOD, fba)
    sbml_setUnits(SMOD, units=None, L3=False)
    C = sbml_setSpeciesL2(SMOD, fba)
    sbml_setCompartmentsL2(SMOD, C)
    sbml_setReactionsL2(SMOD, fba, return_dict=False)
    SMOD.setAnnotation(L3FBA.replace('<?xml version="1.0" encoding="UTF-8"?>', ''))
    sbml_exportSBML2FBAModel(document, fname, directory=directory, return_doc=True)

    if __DEBUG__:
        print('\nModel')
        print(L3FBA)

        import pprint

        PP = pprint.PrettyPrinter()

        print('\nUnits')
        PP.pprint(sbml_setUnits(SMOD, units=None, give_default=True))
        print('\nCompartments')
        c, s = sbml_setSpeciesL2(SMOD, fba, return_dicts=True)
        PP.pprint(c)
        print('\nSpecies')
        PP.pprint(s)
        print('\nReactions')
        PP.pprint(sbml_setReactionsL2(SMOD, fba, return_dict=True))


class FBCconnect(object):
    metaprefix = 'meta_'
    model = None
    doc = None
    sbml = None
    sbmllevel = 3
    sbmlversion = 1
    groupsversion = 1
    sbmlns = None
    fbc = None
    fbcversion = 2
    fbcstrict = True
    maxobjname = ('maximize', 'maximise', 'max')
    minobjname = ('minimize', 'minimise', 'min')
    GROUPS_AVAILABLE = False
    groups = None
    groupList = None

    def __init__(self, fbc_version=2, fbc_strict=True, enable_groups=False):
        """
        Connect SBML Packages

         - *fbc_version* [default=2] fbc version
         - *fbc_strict* [default=True] if using FBC V2 set the strict flag
         - *enable_groups* [default=False] use the Groups package if possible

        """
        assert _HAVE_SBML_, '\nERROR: libSBML required for SBML support'
        self.sbml = libsbml
        self.fbcversion = fbc_version
        assert self._checkPackageRegistry(
            'fbc'
        ), '\nSBML Level 3 FBC package required!\n(http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_%28flux%29)'

        self.sbmlns = self.sbml.SBMLNamespaces(
            self.sbmllevel, self.sbmlversion, 'fbc', self.fbcversion
        )
        self.sbmlns.addNamespace("http://www.w3.org/1999/xhtml", "html")

        if enable_groups and self._checkPackageRegistry('groups'):
            self.sbmlns.addNamespace(
                "http://www.sbml.org/sbml/level3/version1/groups/version1", "groups"
            )
            self.GROUPS_AVAILABLE = True
        else:
            self.GROUPS_AVAILABLE = False

        self.doc = self.sbml.SBMLDocument(self.sbmlns)
        self.doc.setPackageRequired('fbc', False)
        self.model = self.doc.createModel()
        self.fbc = self.model.getPlugin("fbc")

        if self.GROUPS_AVAILABLE:
            s2 = self.doc.setPackageRequired('groups', False)
            self.groups = self.model.getPlugin('groups')
            if s2 != 0 or self.groups == None:
                print('\nWARNING: Groups initialisation error, not enabled.')
                self.GROUPS_AVAILABLE = False

        if fbc_version >= 2 and fbc_strict:
            self.fbcstrict = True
            self.fbc.setStrict(True)

        assert (self.fbc != None), '\nSBML Level 3 FBC package required!\n(http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_%28flux%29)'

    def _checkPackageRegistry(self, pkg):
        pkgs = []
        for i in range(0, self.sbml.SBMLExtensionRegistry.getNumRegisteredPackages()):
            pkgs.append(self.sbml.SBMLExtensionRegistry.getRegisteredPackageName(i))
        return pkg in pkgs

    def createGroup(self, grp):
        """
        Adds a group to the SBML model.

         - *grp* the cbmpy group object

        """
        raise NotImplementedError

    def createGroupLegacy(self, gid, group, name=None):
        """
        Tries to add a group to the model if (assuming the groups package is available).

         - *gid* the group name
         - *groups* a list of strings that refer to model SId's

        """

        gl = self.groups.getListOfGroups()
        if len(gl) > 0:
            gl = [g_.getId() for g_ in gl]
        if gid in gl:
            print('Group {} exists appending items'.format(gid))
            G = self.groups.getGroup(gid)
        else:
            print('Group {} does not exist creating'.format(gid))
            G = self.groups.createGroup()
            G.setId(gid)
            if name != None:
                G.setName(name)
        for m_ in group:
            M = G.createMember()
            out = M.setIdRef(m_)
            if out != 0:
                print(
                    'INFO: cannot add {} to group {} as it is not a valid model id'.format(
                        m_, gid
                    )
                )

        return True

    def createFluxBoundV1(self, fid, reaction, operation, value):
        FB = self.fbc.createFluxBound()
        FB.setId(fid)
        FB.setReaction(reaction)
        FB.setOperation(operation)
        if not numpy.isinf(value):
            FB.setValue(float(value))
        elif value > 0.0:
            FB.setValue(float('inf'))
        elif value < 0.0:
            FB.setValue(float('-inf'))
        else:
            raise RuntimeError(
                "\nCannot interpret bound %s with value %s" % (fid, value)
            )

    #FBC3_VARIABLE_TYPES = {libsbml.FBC_VARIABLE_TYPE_LINEAR : 'linear',
                          #libsbml.FBC_VARIABLE_TYPE_QUADRATIC : 'quadratic',
                          #libsbml.FBC_VARIABLE_TYPE_INVALID : 'invalid',
                          #'linear' : libsbml.FBC_VARIABLE_TYPE_LINEAR,
                          #'quadratic' : libsbml.FBC_VARIABLE_TYPE_QUADRATIC,
                          #'invalid' : libsbml.FBC_VARIABLE_TYPE_INVALID
                          #}


    def createObjective(self, oid, osense, flux_objs, name=None, active=True):
        """
        Create and add the FBC Objective function

         - *oid* objective id
         - *osense* objective sense
         - *flux_objs* [(reaction, coefficient)]
         - *name*
         - *active*

        """
        O = self.fbc.createObjective()
        O.setId(oid)
        if name != None:
            O.setName(name)
        if osense.lower() in self.maxobjname:
            osense = 'maximize'
        elif osense.lower() in self.minobjname:
            osense = 'minimize'
        else:
            raise RuntimeError("\nUnknown objective sense: %s" % osense)
        O.setType(osense)
        if active:
            self.fbc.setActiveObjectiveId(oid)
        for fo_ in flux_objs:
            FO = O.createFluxObjective()
            FO.setReaction(fo_[0])
            FO.setCoefficient(fo_[1])
            if self.fbcversion >= 3:
                FO.setId(fo_[2])
                # assume variable is linear if unknown
                if fo_[3] is None:
                    FO.setVariableType('linear')
                else:
                    print('INFO: creating a type', fo_[3], 'fluxobjective')
                    FO.setVariableType(FBC3_VARIABLE_TYPES[fo_[3]])


    def createGeneAssociationV1(self, rid, assoc, gprid=None):
        """
        Create a gene association for a specified reaction and an association string

         - *rid* a reaction id
         - *assoc* the association string e.g.(b0698 and b0697) or (b0696)

        """
        GPR = self.fbc.createGeneAssociation()
        if gprid == None:
            # try this: get rid of invalid symbols in GPRid

            GPR.setId('{}_gpra'.format(rid))

        else:
            GPR.setId(gprid)

        GA = GPR.createAssociation()
        ass = GA.parseInfixAssociation(assoc)
        if ass == None:
            ret0 = -1
        else:
            ret1 = GPR.setAssociation(ass)
            ret2 = GPR.setReaction(rid)
            if ret1 == 0 and ret2 == 0:
                ret0 = 0
            else:
                ret0 = -1
        if ret0 != 0:
            print(
                'WARNING: Incompatible gene label: \"{}\".({}) only dots and underscores are currently allowed in gene labels'.format(
                    rid, assoc
                )
            )
        else:
            return GPR

    def writeToFile(self, filename):
        self.sbml.writeSBMLToFile(self.doc, str(filename))

    def _cleanUP_(self):
        self.fbc = None
        self.groups = None
        self.model = None
        self.doc = None
        self.sbml = None
        self.sbmlns = None
        self.GROUPS_AVAILABLE = None


# def sbml_writeAnnotations(annotations):
# """
# Writes the annotations dictionary as a BiGG compatible SBML <note> and a new style SBML parameter annotation
# """
# raise NotImplemented
# annoSTR = ''
# for K in annotations:
# val = str(annotations[K]).replace('<','&lt;').replace('>','&gt;')
###  val = str(annotations[K]).replace("'",'').replace('<','&lt;').replace('>','&gt;')
# annoSTR += "<html:p>%s: %s</html:p>\n" % (K, val)
###  print annoSTR
# return annoSTR


class CBMtoSBML3(FBCconnect):
    fba = None
    bound_registry = None
    parameter_map = None
    parameter_cntr = 0

    def __init__(self, fba, fbc_version=2, fbc_strict=True, enable_groups=False):
        """
        Convert a CBM model to SBML level 3 with FBC

        """
        super(CBMtoSBML3, self).__init__(fbc_version, fbc_strict, enable_groups)
        self.fba = fba
        self.parameter_map = {}

        if fba.getId() == '' or fba.getId() == None:
            mid0 = 'CBMPY_Model'
        else:
            mid0 = fba.getId()

        mid = ''
        for l in mid0:
            if l.isalnum():
                mid += l
            else:
                mid += '_'

        if not mid[0].isalpha():
            mid = 'MODEL_' + mid
        self.model.setMetaId(METAPREFIX + mid)
        self.model.setId(mid)

    def addModelAnnotation(self, fba):
        # this should go last
        if fba.miriam != None:
            miriam = fba.miriam.getAllMIRIAMUris()
            if len(miriam) > 0:
                sbml_setCVterms(self.model, miriam, model=True)

    def addModelHistory(self):
        """
        If they exist add CBM model history information to SBML model

        """
        sbmh = self.sbml.ModelHistory()
        GO = False
        if self.fba.DATE_CREATED not in [None, '']:
            sbmh.setCreatedDate(self.sbml.Date(*self.fba.DATE_CREATED))
            GO = True
        if self.fba.DATE_MODIFIED not in [None, '']:
            sbmh.setModifiedDate(self.sbml.Date(*self.fba.DATE_MODIFIED))
            GO = True
        if self.fba.MODEL_CREATORS != None and len(self.fba.MODEL_CREATORS) > 0:
            for c_ in self.fba.MODEL_CREATORS:
                C = self.sbml.ModelCreator()
                if self.fba.MODEL_CREATORS[c_]['firstname'] != None:
                    C.setGivenName(self.fba.MODEL_CREATORS[c_]['firstname'])
                if self.fba.MODEL_CREATORS[c_]['lastname'] != None:
                    C.setFamilyName(self.fba.MODEL_CREATORS[c_]['lastname'])
                if self.fba.MODEL_CREATORS[c_]['organisation'] != None:
                    C.setOrganisation(self.fba.MODEL_CREATORS[c_]['organisation'])
                if self.fba.MODEL_CREATORS[c_]['email'] != None:
                    C.setEmail(self.fba.MODEL_CREATORS[c_]['email'])
                sbmh.addCreator(C)
                del C
            GO = True
        if GO:
            self.model.setModelHistory(sbmh)
        del sbmh

    def addBoundsV2(self, compress_bounds=False, sig_dig=20):
        """
        Add FBC V2 style fluxbounds to model

         - *compress_bounds* [default=False] enable parameter compression
         - *sig_dig* [default=20] round off to significant digits

        """
        if compress_bounds:
            print('INFO: V2 bounds compression enabled')
            bounds = {}
            shared_values = {}
            for r_ in self.fba.reactions:
                print(r_.getLowerBound())
                print(r_.getUpperBound())
                lb = round(r_.getLowerBound(), sig_dig)
                ub = round(r_.getUpperBound(), sig_dig)
                rid = r_.getId()
                bounds[rid] = (lb, ub)

                if lb not in shared_values:
                    shared_values[lb] = ['{}_lb'.format(rid)]
                else:
                    shared_values[lb].append('{}_lb'.format(rid))
                if ub not in shared_values:
                    shared_values[ub] = ['{}_ub'.format(rid)]
                else:
                    shared_values[ub].append('{}_ub'.format(rid))

            bnd_pars = {}
            cntr = 0
            for b in shared_values:
                pid = 'bnd_par{}'.format(cntr)
                self.createFbParameterSharedV2(pid, b)
                bnd_pars[b] = pid
                cntr += 1

            for b_ in bounds:
                if b_ not in self.parameter_map:
                    self.parameter_map[b_] = {}
                    self.parameter_map[b_]['lb'] = bnd_pars[bounds[b_][0]]
                    self.parameter_map[b_]['ub'] = bnd_pars[bounds[b_][1]]
                else:
                    print(
                        'Add parameter has detected duplicate reaction: {} (CBXML:1621)'.format(
                            b_
                        )
                    )
        else:
            for bnd in self.fba.flux_bounds:
                bid = self.createFbParameterV2(bnd, add_cbmpy_anno=True)
                rid = bnd.getReactionId()
                if rid not in self.parameter_map:
                    self.parameter_map[rid] = {}
                btype = bnd.getType()
                if btype == 'lower':
                    self.parameter_map[rid]['lb'] = bid
                elif btype == 'upper':
                    self.parameter_map[rid]['ub'] = bid
                else:
                    self.parameter_map[rid]['lb'] = bid
                    self.parameter_map[rid]['ub'] = bid
                    print('ERROR: strange flux bound assigment type error (CBXML:1635)', btype, self.parameter_map[rid])

    def addBoundsV1(self, autofix=False):
        """
        Add the CBM fluxbounds to SBML

         - *autofix* convert '< >' to '<= >=' default = False

        """
        self.bound_registry = []
        for fb_ in self.fba.flux_bounds:
            rid = None
            if fb_.id not in [None, '']:
                rid = fb_.getId()
            else:
                rid = '%s_%s_bnd' % (fb_.reaction, fb_.is_bound)
            if rid not in self.bound_registry:
                self.bound_registry.append(rid)
                operation = None
                if autofix and fb_.operation == 'less':
                    operation = 'lessEqual'
                elif autofix and fb_.operation == 'greater':
                    operation = 'greaterEqual'
                else:
                    operation = fb_.operation
                    # print 'Illegal operation in bound %s: %s' % (fb_.getId(), fb_.operation)
                self.createFluxBoundV1(rid, fb_.reaction, operation, fb_.value)
            else:
                print('Bound %s already exists, skipping ...' % (rid))

    def addObjectives(self):
        """
        Add the CBM objective function to SBML

        """
        for ob_ in self.fba.objectives:
            active = False
            if ob_.getId() == self.fba.getActiveObjective().getId():
                active = True
            flux_objs = [
                (o2.reaction, float(o2.coefficient), o2.getId(), o2.getType()) for o2 in ob_.flux_objectives
            ]
            self.createObjective(ob_.getId(), ob_.operation, flux_objs, active=active)

    def addGenesV2(
        self,
        parse_from_annotation=False,
        annotation_key='GENE ASSOCIATION',
        add_cbmpy_anno=True,
    ):
        """
        Create genes and add to SBML models (FBC V2)

         - *parse_from_annotation* [default=False] if the gpr_assoc list is empty e.g. the GPR associations have
            not been parsed enabling this will do so.
         - *annotation_key* [default='GENE ASSOCIATION'] the key to use for the gene association if parse_from_annotation is True
           by default it will try one of GENE ASSOCIATION, GENE_ASSOCIATION, gene_association, gene association.
         - *add_cbmpy_anno* [default=True] add PySCeS CBMPy annotation

        """
        if len(self.fba.gpr) == 0 and parse_from_annotation:
            self.fba.createGeneAssociationsFromAnnotations(annotation_key)

        for g in self.fba.genes:
            G = self.fbc.createGeneProduct()

            G.setId(formatSbmlId(g.getId()))
            name = g.getName()
            if name != None:
                G.setName(name)
            G.setMetaId(formatSbmlId(g.getMetaId()))
            G.setLabel(g.getLabel())
            sboterm = g.getSBOterm()
            if sboterm is not None and sboterm != '':
                G.setSBOTerm(str(sboterm))

            notes = g.getNotes()
            if notes != '':
                sbml_setNotes3(G, notes)

            if len(g.annotation) > 0:
                if add_cbmpy_anno:
                    annoSTRnew = sbml_writeKeyValueDataAnnotation(g.annotation)
                    annores = G.appendAnnotation(annoSTRnew)
                    if annores == -3:
                        print('Invalid annotation in Gene', g.getId())
                        print(g.annotation, '\n')
            if g.miriam != None:
                # last blah blah
                sbml_setCVterms(G, g.miriam.getAllMIRIAMUris(), model=False)

    def addGeneProteinAssociationsV1(
        self,
        parse_from_annotation=False,
        annotation_key='GENE ASSOCIATION',
        add_cbmpy_anno=True,
    ):
        """
        Add the CBM geneProtein associations to SBML from the gpr_assoc list (FBC V1)

         - *parse_from_annotation* [default=False] if the gpr_assoc list is empty e.g. the GPR associations have
            not been parsed enabling this will do so.
         - *annotation_key* [default='GENE ASSOCIATION'] the key to use for the gene association if parse_from_annotation is True
           by default it will try one of GENE ASSOCIATION, GENE_ASSOCIATION, gene_association, gene association.
         - *add_cbmpy_anno* add PySCeS CBMPy annotation (still to be implemented)

        """
        if len(self.fba.gpr) == 0 and parse_from_annotation:
            print(
                'addGeneProteinAssociationsV1 is generating gpr\'s from annotation strings'
            )
            self.fba.createGeneAssociationsFromAnnotations(annotation_key)

        for g_ in self.fba.gpr:
            rid = g_.getProtein()
            assoc = g_.getAssociationStr(use_labels=True)
            if rid != None and rid != '':
                if assoc != None and assoc != '':
                    GPR = self.createGeneAssociationV1(rid, assoc, gprid=g_.getId())
                    gprnotes = ''
                    for ge_ in g_.getGeneIds():
                        notes = self.fba.getGene(ge_).getNotes()
                        if notes != '':
                            gprnotes += '{}\n'.format(notes)
                    if gprnotes != '':
                        sbml_setNotes3(GPR, gprnotes)
            else:
                print(
                    'WARNING: Skipping GPR association: \"{}\"\n\"{}\"--> \"{}\"'.format(
                        g_.getId(), rid, assoc
                    )
                )


    def addUserDefinedConstraintsV3(self, add_cbmpy_anno=True):
        """
        Add FBC V3 user defined constraints

        """
        assert (self.fbcversion >= 3), "UserDefinedConstraints support require FBCv3"

        for u in self.fba.user_defined_constraints:
            print(u)
            UC = self.fbc.createUserDefinedConstraint()
            UC.setId(u.getId())
            # TODO - for now lets just assume CBMPy is using values and not parameters
            plb = CBModel.Parameter('udcc_{}_lb'.format(u.getId()), u.getLowerBound())
            pub = CBModel.Parameter('udcc_{}_ub'.format(u.getId()), u.getUpperBound())
            LB = self.createParParameter(plb)
            UB = self.createParParameter(pub)
            self.model.addParameter(LB)
            self.model.addParameter(UB)
            UC.setLowerBound(LB.getId())
            UC.setUpperBound(UB.getId())

            for f in u.constraint_components:
                UCFO = UC.createUserDefinedConstraintComponent()
                UCFO.setId(f.getId())
                # TODO - for now lets just assume CBMPy is using values and not parameters
                #var = CBModel.Parameter('udcc_{}_{}_{}_var'.format(u.getId(), f.getId(), f.getVariable()), f.getCoefficient())
                #VAR = self.createParParameter(var)
                #self.model.addParameter(VAR)
                UCFO.setVariable(f.getVariable())
                UCFO.setCoefficient(f.getCoefficient())
                UCFO.setVariableType(FBC3_VARIABLE_TYPES[f.getType()])



    def createParParameter(self, param, add_cbmpy_anno=True):
        """
        Create a generic SBML parameter from a CBMPy parameter

         - *param* a CBMPy parameter object
         - *add_cbmpy_anno* [default=True] add annotation to SBML object

        """
        # print('createParParameter', param.getId())
        par = self.model.createParameter()
        par.setId(param.getId())
        par.setMetaId('meta_{}'.format(param.getId()))
        par.setName(param.getName())
        par.setValue(param.getValue())
        par.setConstant(True)
        if add_cbmpy_anno:
            if param.getSBOterm() is not None:
                par.setSBOTerm(param.getSBOterm())
            if len(param.annotation) > 0:
                annoSTRnew = sbml_writeKeyValueDataAnnotation(param.annotation)
                annores = par.appendAnnotation(annoSTRnew)
            if param.miriam != None:
                miriam = param.miriam.getAllMIRIAMUris()
                if len(miriam) > 0:
                    sbml_setCVterms(param, miriam, model=False)
        return par

    def createFbParameterV2(self, bnd, add_cbmpy_anno=True):
        """
        Create SBML V2 flux bound parameters for reaction

         - *bnd* object

        """
        # print('createFbParameterV2', bnd.getId())

        par = self.model.createParameter()
        bid = '{}_{}'.format(bnd.getReactionId(), bnd.getType())
        par.setId(bid)
        par.setMetaId('meta_{}'.format(bid))
        name = bnd.getName()
        if name is None:
            name = ''
        par.setName(name)
        par.setValue(bnd.getValue())
        par.setConstant(True)
        par.setSBOTerm('SBO:0000625')

        if add_cbmpy_anno:
            if len(bnd.annotation) > 0:
                annoSTRnew = sbml_writeKeyValueDataAnnotation(bnd.annotation)
                annores = par.appendAnnotation(annoSTRnew)
                if annores == -3:
                    print('Invalid annotation in bound:', bid)
                    # print(bnd.annotation, '\n')
            if bnd.miriam != None:
                miriam = bnd.miriam.getAllMIRIAMUris()
                if len(miriam) > 0:
                    sbml_setCVterms(par, miriam, model=False)
        return bid

    def createFbParameterSharedV2(self, pid, value):
        """
        Create SBML V2 flux bound parameters for reaction

         - *value* parameter value
         - *pid* optional parameter id

         returns parameter id

        """

        par = self.model.createParameter()
        par.setId(pid)
        par.setMetaId('meta_{}'.format(pid))
        name = ''
        par.setName(name)
        par.setValue(value)
        par.setConstant(True)
        par.setSBOTerm('SBO:0000625')


    def _cleanUP_(self):
        self.fba = None
        self.parameter_map = {}
        self.parameter_cntr = 0
        self.fbc = None
        self.groups = None
        self.model = None
        self.doc = None
        self.sbml = None
        self.sbmlns = None


def sbml_writeKeyValueDataAnnotation(annotations):
    """
    Writes the key:value annotations as a KeyValueData annotation (http://pysces.sourceforge.net/KeyValueData)
    """

    annoSTR = (
        '<listOfKeyValueData xmlns="http://pysces.sourceforge.net/KeyValueData">\n'
    )
    for K in annotations:
        Ktype = None
        Kval = None
        if annotations[K] == '' or annotations[K] == None:
            Kval = ""
        else:
            Kval = ESCAPE(str(annotations[K]))

        # Kfix = Kfix.lower()
        Kfix = K

        # annoSTR += ' <data id="%s" type="%s" value="%s"/>\n' % (Kfix, Ktype, Kval)
        annoSTR += ' <data id="{}" value="{}"/>\n'.format(Kfix, Kval)
    annoSTR += '</listOfKeyValueData>\n'
    print(annoSTR)
    return annoSTR


RE_loKVD1 = re.compile('<data.*?/>', re.DOTALL)
RE_loKVD1_id = re.compile('id=".*?"')
# RE_loKVD1_type = re.compile('type=".*?"')
RE_loKVD1_value = re.compile('value=".*?"')


def sbml_readKeyValueDataAnnotation(annotations):
    """
    Reads KeyValueData annotation (http://pysces.sourceforge.net/KeyValueData) and returns a dictionary of key:value pairs
    """

    kvout = {}
    if (
        '<listOfKeyValueData xmlns="http://pysces.sourceforge.net/KeyValueData">'
        in annotations
    ):
        data = re.findall(RE_loKVD1, annotations)
        for D_ in data:
            pid = re.search(RE_loKVD1_id, D_)
            if pid != None:
                pid = pid.group(0).split('=')[1].strip()[1:-1]
            # ptype = re.search(RE_loKVD1_type, D_)
            # if ptype != None:
            # ptype = ptype.group(0).split('=')[1].strip()[1:-1]
            pvalue = re.search(RE_loKVD1_value, D_)
            if pvalue != None:
                pvalue = pvalue.group(0).split('=', 1)[1].strip()[1:-1]
                pvalue = __tagStripper__.unescape(pvalue)
                if pvalue.startswith('[') and pvalue.endswith(']'):
                    try:
                        pvalue = eval(pvalue)
                    except (SyntaxError, NameError):
                        pass

            kvout[pid] = pvalue
    return kvout


def sbml_readFBCv3KeyValuePairs(fbcp):
    """
    Reads FBCv3 KeyValue pair annotation and returns a dictionary of key:value pairs

    - *fbcp* an FBC plugin

    """
    print("sbml_readFBCv3KeyValuePairs coming to a reader near you.", fbcp)

    kv_base = {}
    kv_ext = {}

    for k_ in range(fbcp.getNumKeyValuePairs()):
        kvp = fbcp.getKeyValuePair(k_)
        key = kvp.getKey()
        val =  kvp.getValue()
        kv_base[key] = val
        kv_ext[key] = {'id' : kvp.getId(),
                       'name' : kvp.getName(),
                       'uri' : kvp.getUri(),
        }
        print(kv_base)
        print(kv_ext)

    return kv_base, kv_ext




def sbml_setSpeciesL3(
    model,
    fba,
    return_dicts=False,
    add_cobra_anno=False,
    add_cbmpy_anno=True,
    substance_units=True,
    fbc_version=2
):
    """
    Add the species definitions to the SBML object:

     - *model* and SBML model instance or can be None if *return_dicts* == True
     - *fba* a PySCeSCBM model instance
     - *return_dicts* [default=False] only returns the compartment and species dictionaries without updating the SBML
     - *add_cbmpy_anno* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
     - *add_cobra_anno* [default=False] add COBRA <notes> annotation
     - *substance_units* [default=True] defines the species in amounts rather than concentrations (necessary for default mmol/gdw.h)
     - *fbc_version* [default=2] the FBC version to use

    returns:

     - *compartments* a dictionary of compartments (except when given *return_dicts* argument)

    """

    compartments = fba.getCompartmentIds()

    species = {}
    USE_DEFAULT_COMPARTMENT = False
    DEFAULT_COMPARTMENT = 'Cell'
    for s in fba.species:
        if s.compartment == None or s.compartment == '':
            USE_DEFAULT_COMPARTMENT = True
            s.compartment = DEFAULT_COMPARTMENT
            if DEFAULT_COMPARTMENT not in compartments:
                print(
                    'INFO: Species "{}" has no compartment, creating default "{}".'.format(
                        s.getId(), DEFAULT_COMPARTMENT
                    )
                )
                C = CBModel.Compartment(
                    DEFAULT_COMPARTMENT, DEFAULT_COMPARTMENT, 1.0, 3
                )
                C.setAnnotation('CBMPy_info', 'created by SBML writer')
                s.setAnnotation('CBMPy_info', 'compartment added by SBML writer')
                fba.addCompartment(C)
                compartments = fba.getCompartmentIds()
        elif s.compartment not in compartments:
            print(
                'INFO: Compartment "{}" used by species "{}" is not defined, creating.'.format(
                    s.compartment, s.getId()
                )
            )
            C = CBModel.Compartment(s.compartment, s.compartment, 1.0, 3)
            C.setAnnotation('CBMPy_info', 'created by SBML writer')
            fba.addCompartment(C)
            compartments = fba.getCompartmentIds()
        # print '++', s.charge
        if s.miriam == None:
            miriam = {}
        else:
            miriam = s.miriam.getAllMIRIAMUris()
        species.update(
            {
                s.getId(): {
                    'id': s.getId(),
                    'compartment': s.compartment,
                    'name': s.getName(),
                    'charge': s.charge,
                    'value': s.value,
                    'annotation': s.getAnnotations().copy(),
                    'boundary': s.is_boundary,
                    'chemFormula': s.chemFormula,
                    'miriam': miriam,
                    'sboterm': s.getSBOterm(),
                    'notes': s.getNotes(),
                }
            }
        )

    if return_dicts:
        return compartments, species

    keys = sorted(species)
    # keys.sort()

    for spe in keys:
        s = model.createSpecies()
        s.setId(species[spe]['id'])
        Sfbc = s.getPlugin('fbc')
        # METAID
        s.setMetaId(METAPREFIX + species[spe]['id'])
        s.setName(species[spe]['name'])
        # in theory species are constant at whatever level it is set dX/dT == 0 also there is no way to change them
        s.setConstant(False)  # (2023) Non-boundary conditions species cannot be constant
        # if not (species[spe]['charge'] != None or species[spe]['charge'] != ''):
        # s.setCharge(int(species[spe]['charge']))
        if species[spe]['charge'] not in ['', None, 0,]:  # TODO this needs to be considered - bgoli
            # print species[spe]['charge'], int(species[spe]['charge'])
            if (Sfbc.setCharge(int(species[spe]['charge'])) != libsbml.LIBSBML_OPERATION_SUCCESS):
                print('Unable to set charge for species: {}'.format(species[spe]['id']))
        if species[spe]['chemFormula'] not in ['', None]:
            if (Sfbc.setChemicalFormula(str(species[spe]['chemFormula'])) != libsbml.LIBSBML_OPERATION_SUCCESS):
                print('Unable to set chemFormula for species: {}'.format(species[spe]['id']))

        s.setCompartment(species[spe]['compartment'])

        if species[spe]['boundary']:
            s.setBoundaryCondition(True)
            # s.setConstant(True) # this is technically true but may have other impacts so not enabling
        else:
            s.setBoundaryCondition(False)
            s.setConstant(False)
        # print species[spe]['value'], type(species[spe]['value'])
        if not numpy.isnan(species[spe]['value']) and species[spe]['value'] not in ['', None,]:
            if substance_units:
                # set the species to be in amounts or concentrations default is amounts for mmol/gdw.h
                s.setInitialAmount(float(species[spe]['value']))
            else:
                s.setInitialConcentration(float(species[spe]['value']))
        else: # added 2023 if there are no concentrations defined export the model as inf
            if substance_units:
                # set the species to be in amounts or concentrations default is amounts for mmol/gdw.h
                s.setInitialAmount(numpy.Inf)
            else:
                s.setInitialConcentration(numpy.Inf)
        s.setHasOnlySubstanceUnits(substance_units)

        if len(species[spe]['annotation']) > 0:
            if add_cbmpy_anno:
                if fbc_version < 3:
                    annoSTRnew = sbml_writeKeyValueDataAnnotation(species[spe]['annotation'])
                    annores = s.appendAnnotation(annoSTRnew)
                    if annores == -3:
                        print('Invalid annotation in reaction:', species[spe]['id'])
                        print(species[spe]['annotation'], '\n')
                else:
                    # FBCv3 rules this needs to be extended to deal with new KV pair properties
                    sbml_setFBCv3KeyValuePairs(Sfbc, species[spe]['annotation'])

            if add_cobra_anno:
                annoSTR = sbml_writeAnnotationsAsCOBRANote(
                    species[spe]['annotation']
                )  # GOOD RIDDANCE
                if annoSTR != None:
                    s.setNotes(annoSTR)
            elif species[spe]['notes'] != '' and species[spe]['notes'] is not None:
                sbml_setNotes3(s, species[spe]['notes'])
        if len(species[spe]['miriam']) > 0:
            sbml_setCVterms(s, species[spe]['miriam'], model=False)
        if species[spe]['sboterm'] is not None and species[spe]['sboterm'] != '':
            s.setSBOTerm(str(species[spe]['sboterm']))
        if species[spe]['notes'] is not None and species[spe]['notes'] != '':
            sbml_setNotes3(s, species[spe]['notes'])


def sbml_setReactionsL3Fbc(
    fbcmod, return_dict=False, add_cobra_anno=False, add_cbmpy_anno=True, fbc_version=2):
    """
    Add the FBA instance reactions to the SBML model

     - *fbcmod* a CBM2SBML instance
     - *return_dict* [default=False] if True do not add reactions to SBML document instead return a dictionary description of the reactions
     - *add_cbmpy_anno* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
     - *add_cobra_anno* [default=False] add COBRA <notes> annotation
     - *fbc_version* [default=2] writes either FBC v1 (2013) or v2 (2015) or v3 (2023)

    """

    fba = fbcmod.fba
    gpr_reaction_map = {}
    for gpr in fba.gpr:
        gpr_reaction_map[gpr.getProtein()] = gpr.getId()
    # print(gpr_reaction_map)

    reactions = {}
    for r in fba.reactions:
        reactants = []
        products = []
        for rg in r.getStoichiometry():
            if rg[0] > 0.0:
                products.append(rg)
            else:
                reactants.append(rg)
        if r.miriam == None:
            miriam = {}
        else:
            miriam = r.miriam.getAllMIRIAMUris()
        reactions.update(
            {
                r.getId(): {
                    'id': r.getId(),
                    'reactants': reactants,
                    'products': products,
                    'name': r.getName(),
                    'reversible': r.reversible,
                    'exchange': r.is_exchange,
                    'annotation': r.getAnnotations().copy(),
                    'compartment': r.compartment,
                    'miriam': miriam,
                    'notes': r.getNotes(),
                    'sboterm': r.getSBOterm(),
                    'modifiers': r._modifiers_,
                }
            }
        )
    if return_dict:
        return reactions

    for rxn in reactions:
        # print 'Adding reaction:', reactions[rxn]['id']
        r = fbcmod.model.createReaction()
        if fbc_version >= 2:
            FB = r.getPlugin('fbc')
        r.setId(reactions[rxn]['id'])
        # METAID
        r.setMetaId(METAPREFIX + reactions[rxn]['id'])
        if reactions[rxn]['name'] is not None and reactions[rxn]['name'] != '':
            r.setName(reactions[rxn]['name'])
        r.setFast(False)
        for s in range(len(reactions[rxn]['reactants'])):
            sref = r.createReactant()
            if fbcmod.fbcstrict:
                sref.setConstant(True)
            else:
                sref.setConstant(False)
            sref.setStoichiometry(abs(float(reactions[rxn]['reactants'][s][0])))
            sref.setSpecies(reactions[rxn]['reactants'][s][1])
        for p in range(len(reactions[rxn]['products'])):
            pref = r.createProduct()
            if fbcmod.fbcstrict:
                pref.setConstant(True)
            else:
                pref.setConstant(False)
            pref.setStoichiometry(abs(float(reactions[rxn]['products'][p][0])))
            pref.setSpecies(reactions[rxn]['products'][p][1])

        if reactions[rxn]['reversible']:
            r.setReversible(True)
        else:
            r.setReversible(False)

        if (
            reactions[rxn]['compartment'] is not None
            and reactions[rxn]['compartment'] != ''
        ):
            r.setCompartment(reactions[rxn]['compartment'])

        if len(reactions[rxn]['modifiers']) > 0:
            for mo_ in reactions[rxn]['modifiers']:
                r.addModifier(fbcmod.model.getSpecies(mo_))

        if fbc_version >= 2:
            FB.setLowerFluxBound(fbcmod.parameter_map[reactions[rxn]['id']]['lb'])
            FB.setUpperFluxBound(fbcmod.parameter_map[reactions[rxn]['id']]['ub'])
            if reactions[rxn]['id'] in gpr_reaction_map:
                GPR = fba.getGPRassociation(gpr_reaction_map[reactions[rxn]['id']])
                sbgpr = FB.createGeneProductAssociation()
                sbgpr.setId(GPR.getId())
                if GPR.getName() is not None:
                    sbgpr.setName(GPR.getName())
                # TODO: once I have switched FBCV1 to trees then this can go into operation
                # sbml_createAssociationFromAST(ast.parse(GPR.getAssociationStr()).body[0], sbgpr)
                if GPR.getTree() is not None:
                    sbml_createAssociationFromTreeV2(GPR.getTree(), sbgpr)
                else:
                    pass
                    # print('WARNING: {} cannot create association from tree: {}'.format(GPR.getId(), GPR.getTree()))

                # add annotation
                if len(GPR.annotation) > 0:
                    if add_cbmpy_anno:
                        if fbc_version < 3:
                            annoSTRnew = sbml_writeKeyValueDataAnnotation(GPR.annotation)
                            annores = sbgpr.appendAnnotation(annoSTRnew)
                            if annores == -3:
                                print('Invalid annotation in reaction GPR association', reactions[rxn]['id'])
                        else:
                            # FBCv3 rules this needs to be extended to deal with new KV pair properties
                            sbml_setFBCv3KeyValuePairs(sbgpr, GPR.annotation)

                if GPR.miriam is not None:
                    sbml_setCVterms(sbgpr, GPR.miriam.getAllMIRIAMUris(), model=False)

        # add annotation
        if len(reactions[rxn]['annotation']) > 0:
            if add_cbmpy_anno:
                if fbc_version < 3:
                    annoSTRnew = sbml_writeKeyValueDataAnnotation(reactions[rxn]['annotation'])
                    annores = r.appendAnnotation(annoSTRnew)
                    if annores == -3:
                        print('Invalid annotation in reaction', reactions[rxn]['id'])
                        print(reactions[rxn]['annotation'], '\n')
                else:
                    # FBCv3 rules this needs to be extended to deal with new KV pair properties
                    sbml_setFBCv3KeyValuePairs(FB, reactions[rxn]['annotation'])

        if len(reactions[rxn]['miriam']) > 0:
            sbml_setCVterms(r, reactions[rxn]['miriam'], model=False)

        if reactions[rxn]['sboterm'] is not None and reactions[rxn]['sboterm'] != '':
            r.setSBOTerm(str(reactions[rxn]['sboterm']))

        if add_cobra_anno:
            annoSTR = sbml_writeAnnotationsAsCOBRANote(
                reactions[rxn]['annotation']
            )  # GOOD RIDDANCE
            if reactions[rxn]['notes'] != '' and reactions[rxn]['notes'] is not None:
                annoSTR = '{}\n<br/>\n{}'.format(annoSTR, reactions[rxn]['notes'])
            if annoSTR is not None:
                sbml_setNotes3(r, annoSTR)
        elif reactions[rxn]['notes'] != '' and reactions[rxn]['notes'] is not None:
            sbml_setNotes3(r, reactions[rxn]['notes'])



def sbml_setFBCv3KeyValuePairs(fbcp, kv_pairs):
    """
    Adds key value pairs to any FBCv3 SBase derived object

    - *fbcp* an FBCv3 plugin
    - *kv_pairs* a dictionary of CBMPy Key Value pairs (format will be upgraded over next few versions)

    """
    for kv_ in kv_pairs:
        kvp = fbcp.createKeyValuePair()
        kvp.setKey(kv_)
        if kv_pairs[kv_] is not None:
            kvp.setValue(kv_pairs[kv_])



def sbml_setGroupsL3(cs, fba):
    """
    add groups to the SBML model

     - *cs* a CBMLtoSBML instance
     - *fba* a CBMPy model instance

    """
    if not cs.GROUPS_AVAILABLE:
        print('\nWARNING: Groups package not available.')
        return False

    gids = fba.getGroupIds()
    for grp in fba.groups:
        # global g
        g = cs.groups.createGroup()
        g.setId(str(grp.getId()))
        g.setName(str(grp.getName()))
        metaid = grp.getMetaId()
        # print(grp.getId(), metaid)
        if metaid == None or metaid == '':
            metaid = '{}{}'.format(METAPREFIX, grp.getId())
        g.setMetaId(metaid)

        gkind = grp.getKind()
        if gkind in grp._kinds_:
            if gkind == 'collection':
                g.setKind(cs.sbml.GROUP_KIND_COLLECTION)
            elif gkind == 'partonomy':
                g.setKind(cs.sbml.GROUP_KIND_PARTONOMY)
            elif gkind == 'classification':
                g.setKind(cs.sbml.GROUP_KIND_CLASSIFICATION)
        else:
            g.setKind(cs.sbml.GROUP_KIND_UNKNOWN)

        sbo = grp.getSBOterm()
        if sbo != None:
            g.setSBOTerm(str(sbo))
        notes = grp.getNotes()
        if notes != '' and notes is not None:
            sbml_setNotes3(g, notes)
        lom = g.getListOfMembers()
        metaid = '{}{}_members'.format(METAPREFIX, grp.getId())
        lom.setMetaId(metaid)
        sbo = grp.getSharedSBOterm()
        if sbo != None:
            lom.setSBOTerm(str(sbo))
        notes = grp.getSharedNotes()
        if notes != '' and notes is not None:
            sbml_setNotes3(lom, notes)

        for mid in grp.getMemberIDs():
            m = g.createMember()
            m.setIdRef(str(mid))

        if len(grp.annotation) > 0:
            annoSTRnew = sbml_writeKeyValueDataAnnotation(grp.annotation)
            annores = g.appendAnnotation(annoSTRnew)
            if annores == -3:
                print('Invalid annotation in group:', grp.getId())
                print(grp.annotation, '\n')

        sharedAnno = grp.getSharedAnnotations()
        if len(sharedAnno) > 0:
            annoSTRnew = sbml_writeKeyValueDataAnnotation(sharedAnno)
            annores = lom.appendAnnotation(annoSTRnew)
            if annores == -3:
                print('Invalid annotation in group:', lom.getId())
                print(sharedAnno.annotation, '\n')
        sharedMiriam = grp.getSharedMIRIAMannotations()
        if sharedMiriam != None:
            if len(sharedMiriam) > 0:
                sbml_setCVterms(lom, sharedMiriam, model=False)
        if grp.miriam != None:
            miriam = grp.miriam.getAllMIRIAMUris()
            if len(miriam) > 0:
                sbml_setCVterms(g, miriam, model=False)
    return True


def sbml_writeCOBRASBML(fba, fname, directory=None):
    """
    Takes an FBA model object and writes it to file as a COBRA compatible :

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname

    """
    print(
        '\nWARNING: saving in COBRA format may result in a loss of model information!\n'
    )
    time.sleep(3)
    sbml_writeSBML3FBC(
        fba,
        fname,
        directory,
        autofix=True,
        return_fbc=False,
        sbml_level_version=(3, 1),
        add_groups=False,
    )
    sbml_convertSBML3FBCToCOBRA(fname, fname)
    print('Model exported as: {}'.format(fname))


def sbml_getGeneRefs(association, out):
    """
    Walk through a gene association and extract GeneRefs inspired by Frank

    """
    if isinstance(association, libsbml.GeneProductRef):
        ref = association.getGeneProduct()
        if ref not in out:
            out.append(ref)
    else:
        for i in range(association.getNumAssociations()):
            sbml_getGeneRefs(association.getAssociation(i), out)


def sbml_getGPRasDictFBCv1(node, out):
    """
    Converts a GPR string '((g1 and g2) or g3)' to a gprDict which is returned

     - *node* a Python AST note (e.g. `ast.parse(gprstring).body[0]`)
     - *out* a new dictionary that will be be created in place

    """
    return getGPRasDictFromString(node, out)


def sbml_getGPRasDictFBCv2(association, out, cntr):
    """
    Walk through an SBML L3FBCV2 gene protein association and return a dictionary/tree representation

    """
    key = cntr
    # print(association)
    if isinstance(association, libsbml.GeneProductRef):
        ref = association.getGeneProduct()
        if ref is not None and ref != 'None':
            out[ref] = ref
    elif isinstance(association, libsbml.FbcAnd):
        out['_AND_{}'.format(key)] = {}
        for i in range(association.getNumAssociations()):
            out['_AND_{}'.format(key)].update(
                sbml_getGPRasDictFBCv2(
                    association.getAssociation(i), out['_AND_{}'.format(key)], cntr
                )
            )
            cntr += 1
    elif isinstance(association, libsbml.FbcOr):
        out['_OR_{}'.format(key)] = {}
        for i in range(association.getNumAssociations()):
            out['_OR_{}'.format(key)].update(
                sbml_getGPRasDictFBCv2(
                    association.getAssociation(i), out['_OR_{}'.format(key)], cntr
                )
            )
            cntr += 1
    return out


def sbml_createAssociationFromAST(node, out):
    """
    Converts a GPR string '((g1 and g2) or g3)' to an association via a Python AST.
    In future I will get rid of all the string elements and work only with associations
    and AST's.

     - *node* a Python AST note (e.g. body)
     - *out* a new shiny FBC V2 GeneProductAssociation

    """
    if isinstance(node, ast.Name):
        ref = out.createGeneProductRef()
        ref.setGeneProduct(node.id)

    elif isinstance(node, ast.BinOp):
        left = node.left.id
        right = node.right.id
        ref = out.createGeneProductRef()
        gref = '{}-{}'.format(left, right)
        ref.setGeneProduct(formatSbmlId(gref))
    else:
        if isinstance(node, ast.Expr):
            children = [node.value]
        else:
            children = node.values
        for v in children:
            if isinstance(v, ast.BoolOp) and isinstance(v.op, ast.And):
                newex = out.createAnd()
                # print('And', v)
                # walk(v, newand)
            elif isinstance(v, ast.BoolOp) and isinstance(v.op, ast.Or):
                newex = out.createOr()
                # print('Or', v)
                # walk(v, newor)
            else:
                # print('-->', v)
                newex = out
            sbml_createAssociationFromAST(v, newex)


def sbml_createAssociationFromTreeV2(tree, out):
    """
    Converts a GPR tree to an association

     - *tree* a GPR dict tree
     - *out* a new shiny FBC V2 GeneProductAssociation

    """
    for c in tree:
        if c.startswith('_AND_'):
            newex = out.createAnd()
            sbml_createAssociationFromTreeV2(tree[c], newex)
        elif c.startswith('_OR_'):
            newex = out.createOr()
            sbml_createAssociationFromTreeV2(tree[c], newex)
        else:
            ref = out.createGeneProductRef()
            ref.setGeneProduct(c)
    return out

def sbml_writeSBML3FBC(
    fba,
    fname,
    directory=None,
    sbml_level_version=(3, 1),
    autofix=True,
    return_fbc=False,
    gpr_from_annot=False,
    add_groups=False,
    add_cbmpy_annot=True,
    add_cobra_annot=False,
    xoptions={},
):
    """
    Takes an FBA model object and writes it to file as SBML L3 FBC:

     - *fba* an fba model object
     - *fname* the model will be written as XML to *fname*
     - *directory* [default=None] if defined it is prepended to fname
     - *sbml_level_version* [default=(3,1)] a tuple containing the SBML level and version e.g. (3,1)
     - *autofix* convert <> to <=>=
     - *return_fbc* return the FBC converter instance
     - *gpr_from_annot* [default=False] if enabled will attempt to add the gene protein associations from the annotations
       if no gene protein association objects exist
     - *add_cbmpy_annot* [default=True] add CBMPy KeyValueData annotation. Replaces <notes>
     - *add_cobra_annot* [default=True] add COBRA <notes> annotation
     - *xoptions* extended options

       - *fbc_version* [default=2] write SBML3FBC using version 1 (2013) or version 2 (2015) oe version (2023)
       - *validate* [default=False] validate the output SBML file
       - *compress_bounds* [default=False] try compress output flux bound parameters
       - *zip_model* [default=False] compress the model using ZIP encoding
       - *return_model_string* [default=False] return the SBML XML file as a string


    """
    assert (
        _HAVE_SBML_
    ), "\nSBML not available ... install libSBML with Python bindings for SBML support"

    # load options
    fbc_version = 2
    VALIDATE = False
    compress_bounds = False
    zip_model = False
    return_model_string = False
    if 'fbc_version' in xoptions:
        fbc_version = xoptions['fbc_version']
    if 'validate' in xoptions:
        VALIDATE = xoptions['validate']
    if 'return_model_string' in xoptions:
        return_model_string = xoptions['return_model_string']
    if 'compress_bounds' in xoptions and fbc_version >= 2:
        compress_bounds = xoptions['compress_bounds']
    if 'zip_model' in xoptions:
        zip_model = xoptions['zip_model']
    if fbc_version >= 2:
        autofix = True

    print('\nINFO: using FBC version: {}'.format(fbc_version))

    if fba.getName() in [None, '', ' ']:
        fba.setName('cbmpy_fbc_model')

    fba.setModifiedDate()

    USE_GROUPS = False
    if (
        add_groups
        and hasattr(fba, 'groups')
        and fba.groups != None
        and len(fba.groups) > 0
    ):
        USE_GROUPS = True

    cs3 = CBMtoSBML3(fba, fbc_version, fbc_strict=True, enable_groups=USE_GROUPS)

    if fbc_version == 1:
        cs3.addBoundsV1(autofix=autofix)
    elif fbc_version >= 2:
        cs3.addBoundsV2(compress_bounds=compress_bounds)

    cs3.addObjectives()

    if fbc_version == 1:
        cs3.addGeneProteinAssociationsV1(
            parse_from_annotation=gpr_from_annot,
            annotation_key='GENE_ASSOCIATION',
            add_cbmpy_anno=add_cbmpy_annot,
        )
    elif fbc_version >= 2:
        cs3.addGenesV2(
            parse_from_annotation=gpr_from_annot,
            annotation_key='GENE_ASSOCIATION',
            add_cbmpy_anno=add_cbmpy_annot)

    # create a model
    sbml_setDescription(cs3.model, fba)
    sbml_setUnits(cs3.model, units=None, L3=True)
    sbml_setSpeciesL3(
        cs3.model,
        fba,
        add_cobra_anno=add_cobra_annot,
        add_cbmpy_anno=add_cbmpy_annot,
        substance_units=True,
        fbc_version=fbc_version
    )
    sbml_setCompartmentsL3(cs3.model, fba)
    sbml_setReactionsL3Fbc(
        cs3,
        return_dict=False,
        add_cobra_anno=add_cobra_annot,
        add_cbmpy_anno=add_cbmpy_annot,
        fbc_version=fbc_version
    )
    sbml_setParametersL3Fbc(cs3, add_cbmpy_anno=add_cbmpy_annot, fbc_version=fbc_version)

    if USE_GROUPS:
        sbml_setGroupsL3(cs3, fba)

    # stop libSBML deleting all my annotations
    cs3.addModelHistory()
    cs3.addModelAnnotation(fba)

    if directory != None:
        assert os.path.exists(directory), '\n%s does not exist.' % directory
        fname = os.path.join(directory, fname)

    try:
        UseR = getuser()
    except:
        UseR = ''

    h1 = '<?xml version="1.0" encoding="utf-8"?>\n'
    h1 += (
        '<!-- SBML created with CBMPy ('
        + __version__
        + ') on '
        + time.strftime("%a, %d %b %Y %H:%M:%S")
        + ' -->\n'
    )
    if return_model_string:
        modstr = h1 + cs3.doc.toSBML()
        print('Model returned as string')
    elif zip_model:
        cs3.sbml.writeSBMLToFile(cs3.doc, fname + '.zip')
        print('Model exported as: {}'.format(fname + '.zip'))
    else:
        # raw_input('L%sV%s' % (document.getLevel(),document.getVersion()))
        F = open(fname, 'w')
        F.write(h1 + cs3.doc.toSBML())
        F.flush()
        F.close()
        print('Model exported as: {}'.format(fname))

    # try and keep backwards compatability with FBCv2 hack and FBCv3
    if fba.user_constraints is not None and len(fba.user_constraints) >= 1:
        if fbc_version == 2:
            print("\nWARNING: User defined constraints are included in FBCv3 and you are attempting to save a FBCv2 file. Please save the model using the newer format.\nIn CBMPy 0.9 this will become compulsory.\n")
            time.sleep(5)
            fba.exportUserConstraints(fname + '.user_constraints.json')
        elif fbc_version >= 3:
            print("\nWARNING: User defined constraints defined using old FBCv2 format, converting to FBCv3 ...\n")
            fba.convertUserConstraintsToUserDefinedConstraints()

    # add FBCv3 user defined constraints
    if fbc_version >= 3:
        cs3.addUserDefinedConstraintsV3(add_cbmpy_anno=add_cbmpy_annot)

    if VALIDATE:
        sbml_setValidationOptions(cs3.doc, level='full')
        print('\nPerforming validation on output SBML ...\n')
        errors, warnings, others, DOCUMENT_VALID, MODEL_VALID = sbml_validateDocument(
            cs3.doc, fullmsg=False
        )
        if not DOCUMENT_VALID:
            if return_model_string:
                fname2 = 'model.invalid.xml'
                cs3.sbml.writeSBMLToFile(cs3.doc, fname2)
            else:
                fname2 = fname + '.invalid'
                shutil.move(fname, fname2)
            print('\nSBML document is invalid and written as: {}'.format(fname2))
        if not MODEL_VALID:
            print('\nSBML document is valid but the model is not (or is incomplete):')
            for e in errors:
                print(errors[e])

    else:
        sbml_setValidationOptions(cs3.doc, level='normal')

    cs3._cleanUP_()
    del cs3
    if return_model_string:
        return modstr


def sbml_setValidationOptions(D, level):
    """
    set the validation level of an SBML document

     - *D* an SBML document
     - *level* the level of consistency check can be either one of:

      - 'normal' basic id checking only
      - 'full' all checks enabled
      - None disable all validation

    """

    if level == 'normal':
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_INTERNAL_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, False)
    elif level == 'full':
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_INTERNAL_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, True)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, True)
    elif level is None:
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_INTERNAL_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, False)
        D.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, False)


def sbml_readCOBRASBML(
    fname,
    work_dir=None,
    return_sbml_model=False,
    delete_intermediate=False,
    fake_boundary_species_search=False,
    output_dir=None,
    speciesAnnotationFix=True,
    skip_genes=False,
):
    """
    Read in a COBRA format SBML Level 2 file with FBA annotation where and return either a CBM model object
    or a (cbm_mod, sbml_mod) pair if return_sbml_model=True

     - *fname* is the filename
     - *work_dir* is the working directory
     - *return_sbml_model* [default=False] return a a (cbm_mod, sbml_mod) pair
     - *delete_intermediate* [default=False] delete the intermediate SBML Level 3 FBC file
     - *fake_boundary_species_search* [default=False] after looking for the boundary_condition of a species search for overloaded id's <id>_b
     - *output_dir* [default=None] the directory to output the intermediate SBML L3 files (if generated) default to input directory
     - *speciesAnnotationFix* [default=True]
     - *skip_genes* [default=False] convert GPR associations

    """
    try:
        new_file = sbml_convertCOBRASBMLtoFBC(
            fname, outname=None, work_dir=work_dir, output_dir=output_dir
        )
        print(new_file)
    except Exception as why:
        print('\nCOBRA file conversion failed:\n\"{}\"'.format(why))
        return None
    res = sbml_readSBML3FBC(
        new_file,
        work_dir=work_dir,
        return_sbml_model=False,
        xoptions={'nogenes': skip_genes},
    )
    if fake_boundary_species_search:
        if return_sbml_model:
            cmod = res[0]
        else:
            cmod = res
        for s_ in cmod.species:
            sid = s_.getId()
            if sid[-2:] == '_b':
                print('INFO: Fake boundary (_b) metabolite fixed: {}'.format(sid))
                s_.setBoundary()
    if speciesAnnotationFix:
        if return_sbml_model:
            cmod = res[0]
        else:
            cmod = res
        for s_ in cmod.species:
            try:
                processSpeciesChargeChemFormulaAnnot(
                    s_,
                    getFromName=False,
                    overwriteCharge=True,
                    overwriteChemFormula=True,
                )
            except:
                print(
                    'processSpeciesChargeChemFormulaAnnot failed for species with id: {}'.format(
                        s_.getId()
                    )
                )
        sbml_writeSBML3FBC(cmod, fname=new_file)
    if delete_intermediate:
        os.remove(new_file)
    else:
        print('\nINFO: SBML Level 3 + FBC file generated as: {}'.format(new_file))
    return res


def sbml_convertCOBRASBMLtoFBC(fname, outname=None, work_dir=None, output_dir=None):
    """
    Read in a COBRA SBML Level 2 file and return the name of the created SBML Level 3 with FBC
    file that is created in the output directory

     - *fname* is the filename
     - *outname* the name of the output file. If not specified then <filename>.l3fbc.xml is used as default
     - *work_dir* [default=None] is the working directory
     - *output_dir* [default=None] is the output directory (default is work_dir)

    This method is based on code from libSBML (http://sbml.org) in the file "convertCobra.py"
    written by Frank T. Bergmann.

    """
    assert (
        _HAVE_SBML_
    ), "\nERROR: SBML not available ... install libSBML 5.8.0 r newer with Python bindings for SBML (http://sbml.org/Downloads)"
    assert (
        libsbml.LIBSBML_VERSION >= 50800
    ), "\nWARNING: COBRA SBML import requires libSBML 5.8.0 or newer (your version: {}) (http://sbml.org/Downloads)".format(
        libsbml.LIBSBML_VERSION
    )

    if work_dir != None:
        fname_only = fname
        fname = os.path.join(work_dir, fname)
    else:
        work_dir = os.path.split(fname)
        fname_only = work_dir[1]
        if work_dir[0] != '':
            work_dir = work_dir[0]
        else:
            work_dir = os.getcwd()

    if output_dir == None:
        output_dir = work_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    assert os.path.exists(fname), '\nERROR: file \"{}\" does not exist"'.format(fname)

    if outname != None:
        newfname = outname
    # elif fname[-4:] == '.xml':
    # newfname = fname_only.replace('.xml', '.l3fbc.xml')
    else:
        newfname = fname_only + '.l3fbc.xml'

    newfname = os.path.join(output_dir, newfname)

    print('\nwork_dir: {}'.format(work_dir))
    print('output_dir: {}'.format(output_dir))
    print('fname: {}'.format(fname))
    print('newfname: {}'.format(newfname))

    print('\nRead ...')
    sbmldoc = SBMLreader.readSBML(fname)
    print('Read reports {} errors'.format(sbmldoc.getNumErrors()))

    if sbmldoc.getNumErrors() > 0:
        if sbmldoc.getError(0).getErrorId() == libsbml.XMLFileUnreadable:
            # Handle case of unreadable file here.
            sbmldoc.printErrors()
        elif sbmldoc.getError(0).getErrorId() == libsbml.XMLFileOperationError:
            # Handle case of other file error here.
            sbmldoc.printErrors()
        else:
            # Handle other error cases here.
            sbmldoc.printErrors()
        return None

    props = libsbml.ConversionProperties()
    props.addOption("convert cobra", True, "convert cobra sbml to fbc")
    print('\nConvert ...')
    result = sbmldoc.convert(props)
    print('Convert returns result {}'.format(result))
    if result != libsbml.LIBSBML_OPERATION_SUCCESS:
        print("[Error] Conversion failed... (%d)" % (result))
        return None
    SBMLwriter.writeSBML(sbmldoc, str(newfname))
    print("\nINFO: successfully converted file {} to {}\n".format(fname, newfname))
    props = None
    del sbmldoc
    sbmldoc = None
    return newfname


def sbml_convertSBML3FBCToCOBRA(fname, outname=None, work_dir=None, output_dir=None):
    """
    Read in a SBML Level 3 file and return the name of the created COBRA
    file that is created in the output directory

     - *fname* is the filename
     - *outname* the name of the output file. If not specified then <filename>.cobra.xml is used as default
     - *work_dir* [default=None] is the working directory
     - *output_dir* [default=None] is the output directory (default is work_dir)

    This method is based on code from libSBML (http://sbml.org) in the file "convertFbcToCobra.py"
    written by Frank T. Bergmann.

    """
    assert (
        _HAVE_SBML_
    ), "\nERROR: SBML not available ... install libSBML 5.8.0 r newer with Python bindings for SBML (http://sbml.org/Downloads)"
    assert (
        libsbml.LIBSBML_VERSION >= 50800
    ), "\nWARNING: COBRA SBML import requires libSBML 5.8.0 or newer (your version: {}) (http://sbml.org/Downloads)".format(
        libsbml.LIBSBML_VERSION
    )

    # global SBMLreader, SBMLwriter

    if work_dir != None:
        fname_only = fname
        fname = os.path.join(work_dir, fname)
    else:
        work_dir = os.path.split(fname)
        fname_only = work_dir[1]
        if work_dir[0] != '':
            work_dir = work_dir[0]
        else:
            work_dir = os.getcwd()

    if output_dir == None:
        output_dir = work_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    assert os.path.exists(fname), '\nERROR: file \"{}\" does not exist"'.format(fname)

    if outname != None:
        newfname = outname
    # elif fname[-4:] == '.xml':
    # newfname = fname_only.replace('.xml', '.l3fbc.xml')
    else:
        newfname = fname_only + '.cobra.xml'

    newfname = os.path.join(output_dir, newfname)

    print('\nwork_dir: {}'.format(work_dir))
    print('output_dir: {}'.format(output_dir))
    print('fname: {}'.format(fname))
    print('newfname: {}'.format(newfname))

    # SBMLreader  = libsbml.SBMLReader()

    print('\nRead ...')
    sbmldoc = SBMLreader.readSBML(fname)
    print(
        'Converter reports {} errors, this is probably normal.'.format(
            sbmldoc.getNumErrors()
        )
    )
    # debug stuff
    # M = sbmldoc.getModel()
    # for s_ in range(M.getNumSpecies()):
    # s = M.getSpecies(s_)
    # s.unsetAnnotation()
    # s.unsetNotes()
    # for r_ in range(M.getNumReactions()):
    # r = M.getReaction(r_)
    # r.unsetAnnotation()
    # r.unsetNotes()

    # if sbmldoc.getNumErrors() > 0:
    # if sbmldoc.getError(0).getErrorId() == libsbml.XMLFileUnreadable:
    ## Handle case of unreadable file here.
    # sbmldoc.printErrors()
    # elif sbmldoc.getError(0).getErrorId() == libsbml.XMLFileOperationError:
    ## Handle case of other file error here.
    # sbmldoc.printErrors()
    # else:
    ## Handle other error cases here.
    # sbmldoc.printErrors()
    # return None

    props = libsbml.ConversionProperties()
    props.addOption("convert fbc to cobra", True, "Convert FBC model to Cobra model")
    print('\nConvert ...')
    result = sbmldoc.convert(props)
    print('Convert returns result {}'.format(result))
    if result != libsbml.LIBSBML_OPERATION_SUCCESS:
        print("[Error] Conversion failed... (%d)" % (result))
        return None
    SBMLwriter.writeSBML(sbmldoc, str(newfname))
    print("\nINFO: successfully converted file {} to {}\n".format(fname, newfname))
    del sbmldoc
    return newfname


def sbml_validateDocument(D, fullmsg=False, docread=False):
    """
    Validates and SBML document returns three dictionaries, errors, warnings, other and a boolean indicating an invalid document:

     - *D* and SBML document
     - *fullmsg* [default=False] optionally display the full error message

    """

    warnings = {}
    errors = {}
    others = {}
    DOCUMENT_VALID = True
    MODEL_VALID = True
    e = None
    if docread and D.getNumErrors() > 0:
        sbmlerrors = D.getNumErrors()
    else:
        sbmlerrors = D.validateSBML()

    for e_ in range(sbmlerrors):
        e = D.getError(e_)
        eid = e.getErrorId()
        msg = {
            'line': [e.getLine()],
            'msg': e.getShortMessage(),
            'fullmsg': e.getMessage(),
            'severity_int': e.getSeverity(),
            'package': e.getPackage(),
            'severity': e.getSeverityAsString(),
            'isfatal': e.isFatal(),
            'isxml': e.isXML(),
        }
        if msg['isfatal']:
            DOCUMENT_VALID = False
            MODEL_VALID = False
        # fail model validaty on model error
        if msg['severity_int'] == 2:
            MODEL_VALID = False

        if e.isError():
            if eid not in errors:
                errors[eid] = msg
            else:
                if msg['line'][0] not in errors[eid]['line']:
                    errors[eid]['line'].append(msg['line'][0])
        elif e.isWarning():
            if eid not in warnings:
                warnings[eid] = msg
            else:
                if msg['line'][0] not in warnings[eid]['line']:
                    warnings[eid]['line'].append(msg['line'][0])

        else:
            if eid not in others:
                others[eid] = msg
            else:
                if msg['line'][0] not in others[eid]['line']:
                    others[eid]['line'].append(msg['line'][0])
    del e

    not_relevant = [20221, 20616, 80601, 99130, 99508]
    for er in warnings:
        if er in not_relevant:
            warnings[er]['line'] = []
            warnings[er]['severity'] = 'NAtoFBC'
        elif warnings[er]['line'] == [1]:
            warnings[er]['line'] = []
        if fullmsg:
            warnings[er]['msg'] = warnings[er]['fullmsg']
            warnings[er]['fullmsg'] = ''

    print('Validation report:\n==================\n')
    if len(errors) > 0:
        print('\nErrors\n------\n')
        eidx = list(errors.keys())
        eidx.sort()
        for e_ in eidx:
            print(
                'Error {} (severity={}):\n\n - {}\n - lines: {}\n'.format(
                    e_, errors[e_]['severity'], errors[e_]['msg'], errors[e_]['line']
                )
            )
    if len(warnings) > 0:
        print('\nWarnings\n--------\n')
        eidx = list(warnings.keys())
        eidx.sort()
        for e_ in eidx:
            print(
                'Warning {} (severity={}):\n\n - {}\n - lines: {}\n'.format(
                    e_,
                    warnings[e_]['severity'],
                    warnings[e_]['msg'],
                    warnings[e_]['line'],
                )
            )
    if len(others) > 0:
        print('\nOther\n-----\n')
        eidx = list(others.keys())
        eidx.sort()
        for e_ in eidx:
            print(
                'Info {} (severity={}):\n\n - {}\n - lines: {}\n'.format(
                    e_, others[e_]['severity'], others[e_]['msg'], others[e_]['line']
                )
            )
    print('End.\n')

    return errors, warnings, others, DOCUMENT_VALID, MODEL_VALID


def sbml_fileValidate(f, level="normal"):
    """
    Validate an SBML file and model

     - *f* the SBML file
     - *level* [default='normal'] the level of validation "normal" or "full"

    """
    if not os.path.exists(f):
        print('ERROR: invalid file')
        return None, None
    D = libsbml.readSBMLFromFile(f)
    sbml_setValidationOptions(D, level=level)
    print('\nPerforming validation on input SBML ...\n')
    DOCUMENT_VALID = True
    errors = warnings = others = {}
    errors, warnings, others, DOCUMENT_VALID, MODEL_VALID = sbml_validateDocument(D)
    del D
    if not DOCUMENT_VALID:
        print('\nValidation \"{}\" has detected an invalid SBML document'.format(level))
    elif not MODEL_VALID:
        print('\nValidation \"{}\" has detected an invalid SBML model'.format(level))
    else:
        print('\nValidation \"{}\" successful'.format(level))
    return DOCUMENT_VALID, MODEL_VALID, errors


def sbml_fileFindVersion(f):
    """
    Try and find the SBML version and FBC support

     - *f* the SBML file

    """
    if not os.path.exists(f):
        print('ERROR: invalid file')
        return None
    D = libsbml.readSBMLFromFile(f)
    ns = D.getNamespaces()
    uris = []
    for n in range(ns.getNumNamespaces()):
        uris.append(ns.getURI(n))
    output = None
    for ns, idx in SBML_NS:
        if ns in uris:
            output = idx
            break
    if output == 'L2':
        F = open(f, 'r')
        l2type = None
        for l in F:
            if (
                '<fba:fluxBalance xmlns:fba="http://www.sbml.org/sbml/level3/version1/fba/version1">'
                in l
            ):
                l2type = 'L2FBA'
                break
            if '<parameter id="OBJECTIVE_COEFFICIENT" value="0"' in l:
                l2type = 'COBRA'
                break
        F.close()
        if l2type is not None:
            output = l2type

    msg = ''
    if output == 'L3V1FBC1':
        msg = 'SBML Level 3 FBC version 1 model detected, loading with cbmpy.readSBML3FBC()'
    elif output == 'L3V1FBC2':
        msg = 'SBML Level 3 FBC version 2 model detected, loading with cbmpy.readSBML3FBC()'
    elif output == 'L3V1FBC3':
        msg = 'Awesome! SBML Level 3 FBC version 3 model detected, loading with cbmpy.readSBML3FBC()'
    elif output == 'L2FBA':
        msg = 'SBML Level 2 FAME model detected, loading with cbmpy.readSBML2FBA()'
    elif output == 'COBRA':
        msg = 'COBRA SBML L2 model detected, loading with cbmpy.readCOBRASBML()'
    else:
        msg = 'Unknown model type, contact the SBML community or developer for more details.'
    del D
    print(msg)
    return output, msg


def setCBSBOterm(sbo, obj):
    """
    Given an SBOterm from libSBML, add it to a CBMPy object

     - *sbo* the sbo term string
     - *obj* the CBMPy Fbase derived object

    """
    if sbo is not None and sbo != '':
        try:
            obj.setSBOterm(sbo)
        except AssertionError:
            print(
                'WARNING: {} is not a valid SBO term and is being ignored.'.format(sbo)
            )


# class decorator, from six, that will be used to solve Python2/3 metaclass issues
# usage is NewClass(with_metaclass(MetaClass, BaseClasses)):


def with_metaclass(meta, *bases):
    """
    Create a base class with a metaclass.
    Usage is: NewClass(with_metaclass(MetaClass, BaseClass*)
    """

    # This requires a bit of explanation: the basic idea is to make a dummy
    # metaclass for one level of class instantiation that replaces itself with
    # the actual metaclass.

    class metaclass(type):
        def __new__(cls, name, this_bases, d):
            return meta(name, bases, d)

        @classmethod
        def __prepare__(cls, name, this_bases):
            return meta.__prepare__(name, bases)

    return type.__new__(metaclass, 'temporary_class', (), {})


def sbml_readSBML3FBC(fname, work_dir=None, return_sbml_model=False, xoptions={}):
    """
    Read in an SBML Level 3 file with FBC annotation where and return either a CBM model object
    or a (cbm_mod, sbml_mod) pair if return_sbml_model=True

     - *fname* is the filename
     - *work_dir* is the working directory
     - *return_sbml_model* [default=False] return a a (cbm_mod, sbml_mod) pair
     - *xoptions* special load options enable with option = True

       - *nogenes* do not load/process genes
       - *noannot* do not load/process any annotations
       - *validate* validate model and display errors and warnings before loading
       - *readcobra* read the cobra annotation
       - *read_model_string* [default=False] read the model from a string (instead of a filename) containing an SBML document
       - *nmatrix_type* [default='normal'] define the type of stoichiometrich matrix to be built
       - *model_extension_class* extend CBModel class with new class (experimental, Python 3 only)
       - *model_metaclass* add a custom metaclass to CBModel (experimental, Python 3 only)

         - 'numpy' dense numpy array (best performance)
         - 'scipy_csr' scipy sparse matrix (lower performance, low memory)
         - 'sympy' a sympy rational matrix (low performance, high memory, cast to dense to analyse)
         - None do not build matrix

    """

    time00 = time.time()
    assert (
        _HAVE_SBML_
    ), "\nSBML not available ... install libSBML with Python bindings for SBML support"

    # load options
    LOADGENES = True
    LOADANNOT = True
    DEBUG = False
    VALIDATE = False
    READCOBRA = False
    READ_MODEL_STRING = False
    NMATRIX_TYPE = 'numpy'
    CUSTOM_MODEL_EXTENSION = False
    CUSTOM_MODEL_METACLASS = False

    if 'nogenes' in xoptions and xoptions['nogenes']:
        LOADGENES = False
        print('\nGPR loading disabled!\n')
    if 'noannot' in xoptions and xoptions['noannot']:
        LOADANNOT = False
        print('\nAnnotation loading disabled!\n')
    if 'debug' in xoptions and xoptions['debug']:
        DEBUG = True
        print('\nDebug enabled!\n')
    if 'validate' in xoptions and xoptions['validate']:
        VALIDATE = True
    if 'readcobra' in xoptions and xoptions['readcobra']:
        READCOBRA = True
    if 'read_model_string' in xoptions and xoptions['read_model_string']:
        READ_MODEL_STRING = True
    if 'nmatrix_type' in xoptions:
        if xoptions['nmatrix_type'] is None or xoptions['nmatrix_type'] == 'None':
            NMATRIX_TYPE = None
        elif xoptions['nmatrix_type'] == 'numpy':
            NMATRIX_TYPE = 'numpy'
        elif xoptions['nmatrix_type'] == 'scipy_csr':
            NMATRIX_TYPE = 'scipy_csr'
        elif xoptions['nmatrix_type'] == 'sympy':
            NMATRIX_TYPE = 'sympy'
    if 'model_extension_class' in xoptions:
        if os.sys.version_info > (3, 0):
            CUSTOM_MODEL_EXTENSION = True

            class CBModelExtended(CBModel.Model, xoptions['model_extension_class']):
                __CUSTOM_MODEL_EXTENSION__ = xoptions['model_extension_class'].__name__

                def __init__(self, pid):
                    super().__init__(pid)

        else:
            print(
                '\nWARNING: Custom CBModel extension is supported in Python 3 only.\n'
            )

    if 'model_metaclass' in xoptions:
        raise NotImplementedError("Custom metaclass loading is still being considered.")
        # if os.sys.version_info > (3, 0):
        # CUSTOM_MODEL_METACLASS = True
        # class CBModelExtended(CBModel.Model, metaclass=xoptions['model_metaclass']):
        # __CUSTOM_MODEL_METACLASS__ = xoptions['model_metaclass'].__name__
        # def __init__(self, pid):
        # super().__init__(pid)
        # else:
        # print('\nWARNING: Custom CBModel metaclass is supported in Python 3 only.\n')

    D = None
    if READ_MODEL_STRING:
        D = libsbml.readSBMLFromString(fname)
        fname = 'string_source'
    else:
        if work_dir != None:
            D = libsbml.readSBMLFromFile(str(os.path.join(work_dir, fname)))
        else:
            D = libsbml.readSBMLFromFile(str(fname))

    # check for file read errors
    if D.getNumErrors() > 0:
        print('***\nDANGER: libSBML reports *{}* read errors, this indicates that this is an invalid SBML file. I will try to load it but cannot guarantee its accuracy.\n***\n'.format(
                D.getNumErrors()))
        errors, warnings, others, DOCUMENT_VALID, MODEL_VALID = sbml_validateDocument(D, docread=True)
        if not DOCUMENT_VALID:
            raise RuntimeError(
                "\nERROR: Validation has detected an invalid SBML document This is a fatal error.\n"
            )
        if not MODEL_VALID:
            print(
                '\n***\nERROR: SBML model contains errors but I will try to load it anyway. Depending on the severity of the error I cannot guarantee it\'s accuracy or consistency!!.\n***\n'
            )
        if not DOCUMENT_VALID or not MODEL_VALID:
            time.sleep(5)
    # set consistency checking level for document

    if VALIDATE:
        sbml_setValidationOptions(D, level='full')
        print('\nPerforming validation on input SBML ...\n')
        errors, warnings, others, DOCUMENT_VALID, MODEL_VALID = sbml_validateDocument(D)
        if not DOCUMENT_VALID:
            raise RuntimeError("\nValidation has detected an invalid SBML document")
        elif not MODEL_VALID:
            for e in errors:
                print(errors[e])
            print('\nERROR: SBML document is valid but the model may contain errors\n')
            time.sleep(2)
    else:
        sbml_setValidationOptions(D, level='normal')

    M = D.getModel()
    assert M != None, "\n\nInvalid SBML file"
    assert (M.getLevel() >= 3 and M.getVersion() >= 1), "\nAn SBML L3V1 or greater model is required"

    ## we are now going to allow loading non-FBC models (just to be friendly)
    HAVE_FBC = True
    if not hasattr(libsbml, 'Objective'):
        print(
            "\n\nSBML FBC package required, see http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_%28flux%29 for more details."
        )
        HAVE_FBC = False
    if HAVE_FBC:
        FBCplg = M.getPlugin(str("fbc"))
        if FBCplg is None:
            HAVE_FBC = False
    if not HAVE_FBC:
        print(
            '\nModel is not an SBML3 FBC model.\nPlease try cbmpy.readCOBRASBML(\'{}\') for models encoded in the COBRA dialect\nor cbmpy.readSBML2FBA(\'{}\') for models in FAME format.\n'.format(
                fname, fname
            )
        )
        # return None

    fbc_version = 0
    FBCstrict = True
    if HAVE_FBC:
        fbc_version = FBCplg.getPackageVersion()
        if fbc_version >= 2:
            FBCstrict = FBCplg.getStrict()

        if fbc_version >= 2 and not FBCstrict:
            print("\nWARNING!!!!\n")
            print("This model has fbc:strict=\"false\" this means that is not necessarily a linear program and may contain a number of unsupported features containing aribtrary mathematical expressions such as, InitialAssignments, Rules, Events etc.")
            print("\nCBMPy can continue to load this model but will treat it as a convex linear problem and only load what it can interpret.")
            print("\nWARNING!!!!\n")
            if not input('\nDo you wish to continue (Y/N): ') == 'Y':
                os.sys.exit(-1)

    # print some model information
    print('FBC version: {}'.format(fbc_version))
    print('M.getNumReactions: {}'.format(M.getNumReactions()))
    print('M.getNumSpecies: {}'.format(M.getNumSpecies()))
    if HAVE_FBC:
        print('FBC.getNumObjectives: {}'.format(FBCplg.getNumObjectives()))
        if fbc_version == 1:
            print('FBC.getNumGeneAssociations: {}'.format(FBCplg.getNumGeneAssociations()))
            print('FBC.getNumFluxBounds: {}'.format(FBCplg.getNumFluxBounds()))
        elif fbc_version >= 2:
            print('FBC.getNumParameters: {}'.format(M.getNumParameters()))
            print('FBC.getNumGeneProducts: {}'.format(FBCplg.getNumGeneProducts()))

    model_id = M.getId()
    if model_id in [None, '']:
        model_id = 'sbml_model'
    model_name = M.getName()
    model_description = sbml_getNotes(M)
    # this has been labelled as too agressive and temporarily removed
    # model_description = xml_stripTags(model_description).strip()
    model_description = model_description.strip()

    # TODO bgoli this will need to be refactored and scrubbed at some point
    __HAVE_FBA_ANOT_OBJ__ = True
    __HAVE_FBA_ANOT_BNDS__ = True
    __HAVE_FBA_ANOT_GENEASS__ = True
    if HAVE_FBC:
        if FBCplg.getNumFluxBounds() < 1:
            __HAVE_FBA_ANOT_BNDS__ = False
        if FBCplg.getNumObjectives() < 1:
            __HAVE_FBA_ANOT_OBJ__ = False
        if FBCplg.getNumGeneAssociations() < 1:
            __HAVE_FBA_ANOT_GENEASS__ = False
    else:
        __HAVE_FBA_ANOT_OBJ__ = False
        __HAVE_FBA_ANOT_BNDS__ = False
        __HAVE_FBA_ANOT_GENEASS__ = False

    time0 = time.time()

    SPEC = []
    for s in range(M.getNumSpecies()):
        # add support for initialAmount and hasOnlySubstanceUnits
        SBSp = M.getSpecies(s)
        boundCon = False
        if SBSp.getBoundaryCondition():
            if __DEBUG__:
                print('Real boundary metabolite: {}'.format(SBSp.getId()))
            boundCon = True
        CF = None  # chemical formula
        CH = None

        if HAVE_FBC:
            SBSpF = SBSp.getPlugin("fbc")
            if SBSpF != None:
                CF = SBSpF.getChemicalFormula()
                if fbc_version < 3:
                    CH = int(SBSpF.getCharge())
                else:
                    CH = float(SBSpF.getCharge())
            # print CF, CH

        NM = SBSp.getName()  # get name
        # to strip a BiGG file see CBTools
        S = CBModel.Species(
            SBSp.getId(),
            boundary=boundCon,
            name=NM,
            value=SBSp.getInitialConcentration(),
            compartment=SBSp.getCompartment(),
            charge=CH,
            chemFormula=CF,
        )
        S.annotation = {}
        if LOADANNOT:
            if fbc_version < 3:
                S.annotation = sbml_readKeyValueDataAnnotation(SBSp.getAnnotationString())
                if S.annotation == {}:
                    if libsbml.getLibSBMLVersion() >= 51903:
                        node_txt = libsbml.XMLNode.convertXMLNodeToString(SBSp.getNotes())
                    else:
                        node_txt = libsbml.XMLNode_convertXMLNodeToString(SBSp.getNotes())
                    S.annotation = sbml_readCOBRANote(node_txt)
                    # SBSp.unsetNotes()
            else:
                # TODO bgoli deal with v3 extended annotation
                S.annotation, S.annotation_ext = sbml_readFBCv3KeyValuePairs(SBSpF)

            manot = sbml_getCVterms(SBSp, model=False)
            if manot != None:
                S.miriam = manot
            del manot

        setCBSBOterm(SBSp.getSBOTermID(), S)
        S.setNotes(sbml_getNotes(SBSp))
        SPEC.append(S)

    boundary_species = [s.getId() for s in SPEC if s.is_boundary]
    spec_id = [s.getId() for s in SPEC]

    if DEBUG:
        print('Species load: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # we'll now deal with global parameters even if they are not used
    PARAM_D = {}
    for p_ in range(M.getNumParameters()):
        P = M.getParameter(p_)
        pid = P.getId()
        pdict = {
            'id': pid,
            'value': P.getValue(),
            'constant': P.getConstant(),
            'sbo': P.getSBOTermID(),
            'name': P.getName(),
            'annotation': None,
            'annotation_ext': None,
            'miriam': None,
            'association': [],
            'notes': sbml_getNotes(P),
            'is_fluxbound': False,
        }
        if LOADANNOT:
            if fbc_version < 3:
                pdict['annotation'] = sbml_readKeyValueDataAnnotation(P.getAnnotationString())
            else:
                # TODO bgoli deal with v3 extended annotation
                Pfbc = P.getPlugin("fbc")
                pdict['annotation'], pdict['annotation_ext'] = sbml_readFBCv3KeyValuePairs(Pfbc)

            manot = sbml_getCVterms(P, model=False)

            if manot != None:
                pdict['miriam'] = manot
            del manot
        PARAM_D[pid] = pdict

    GENE_D = {}
    if HAVE_FBC and fbc_version >= 2:
        for g_ in range(FBCplg.getNumGeneProducts()):
            G = FBCplg.getGeneProduct(g_)
            gid = G.getId()
            gdict = {
                'id': gid,
                #'value' : G.getValue(),
                #'constant' : P.getConstant(),
                'sbo': G.getSBOTermID(),
                'name': G.getName(),
                'label': G.getLabel(),
                'annotation': None,
                'miriam': None,
                'notes': sbml_getNotes(G),
            }
            if LOADANNOT:
                if fbc_version < 3:
                    gdict['annotation'] = sbml_readKeyValueDataAnnotation(G.getAnnotationString())
                else:
                    # TODO bgoli deal with v3 extended annotation
                    Gfbc = G.getPlugin("fbc")
                    gdict['annotation'], gdict['annotation_ext'] = sbml_readFBCv3KeyValuePairs(Gfbc)

                manot = sbml_getCVterms(G, model=False)

                if manot != None:
                    gdict['miriam'] = manot
                del manot
            GENE_D[gid] = gdict

    REAC = []
    FB_data = []
    GPR_D = {}
    reactionIDs = []
    reactionsReversability = []
    gprregex = re.compile('\w+')

    for r in range(M.getNumReactions()):
        SBRe = M.getReaction(r)
        R_id = SBRe.getId()
        if HAVE_FBC and fbc_version >= 2:
            # deal with new style fluxbounds
            RFBCplg = SBRe.getPlugin('fbc')
            lfbid = RFBCplg.getLowerFluxBound()
            ufbid = RFBCplg.getUpperFluxBound()
            if lfbid != '':
                fbl = {
                    'reaction': R_id,
                    'operation': 'greaterEqual',
                    'value': PARAM_D[lfbid]['value'],
                    'id': '{}_lb'.format(R_id),
                    'parameter': lfbid,
                    'annotation': PARAM_D[lfbid]['annotation'],
                    'annotation_ext': PARAM_D[lfbid]['annotation_ext'],
                    'miriam': PARAM_D[lfbid]['miriam'],
                    'sbo': PARAM_D[lfbid]['sbo'],
                    'type': 'lower',
                    'name': PARAM_D[lfbid]['name'],
                }
                PARAM_D[lfbid]['association'].append(R_id)
                PARAM_D[lfbid]['is_fluxbound'] = True
                FB_data.append(fbl)
            if ufbid != '':
                fbu = {
                    'reaction': R_id,
                    'operation': 'lessEqual',
                    'value': PARAM_D[ufbid]['value'],
                    'id': '{}_ub'.format(R_id),
                    'parameter': ufbid,
                    'annotation': PARAM_D[ufbid]['annotation'],
                    'annotation_ext': PARAM_D[ufbid]['annotation_ext'],
                    'miriam': PARAM_D[ufbid]['miriam'],
                    'sbo': PARAM_D[ufbid]['sbo'],
                    'type': 'upper',
                    'name': PARAM_D[ufbid]['name'],
                }
                PARAM_D[ufbid]['association'].append(R_id)
                PARAM_D[ufbid]['is_fluxbound'] = True
                FB_data.append(fbu)

            # deal with new gene associations (why larry why ...)
            SBgpr = RFBCplg.getGeneProductAssociation()
            if SBgpr != None:
                GPR_id = SBgpr.getId()
                if GPR_id == '' or GPR_id is None:
                    GPR_id = '{}_gpr'.format(R_id)
                ass = SBgpr.getAssociation()
                if ass != None:
                    GPR_D[GPR_id] = {'gpr_by_id': ass.toInfix()}
                else:
                    GPR_D[GPR_id] = {'gpr_by_id': ''}
                GPR_D[GPR_id]['reaction'] = R_id

                GPR_D[GPR_id]['gene_ids'] = []
                ## the smart way
                if SBgpr.getAssociation() is not None:
                    sbml_getGeneRefs(SBgpr.getAssociation(), GPR_D[GPR_id]['gene_ids'])
                    GPR_D[GPR_id]['gpr_tree'] = sbml_getGPRasDictFBCv2(
                        SBgpr.getAssociation(), {}, 0
                    )

                gene_ids_sorted = sorted(GPR_D[GPR_id]['gene_ids'], key=len)
                gene_ids_sorted.reverse()
                GPR_D[GPR_id]['gpr_by_label'] = GPR_D[GPR_id]['gpr_by_id']
                GPR_D[GPR_id]['gpr_by_name'] = GPR_D[GPR_id]['gpr_by_id']
                GPR_D[GPR_id]['gene_labels'] = []
                gene_labels_sorted = [GENE_D[x_]['label'] for x_ in gene_ids_sorted]
                gene_names_sorted = [GENE_D[x_]['name'] for x_ in gene_ids_sorted]
                for x_ in range(len(gene_labels_sorted)):
                    if libsbml.LIBSBML_VERSION >= 51300:
                        GPR_D[GPR_id]['gpr_by_id'] = GPR_D[GPR_id]['gpr_by_id'].replace(
                            gene_labels_sorted[x_], gene_ids_sorted[x_]
                        )
                        GPR_D[GPR_id]['gpr_by_name'] = GPR_D[GPR_id][
                            'gpr_by_name'
                        ].replace(gene_labels_sorted[x_], gene_names_sorted[x_])
                    else:
                        GPR_D[GPR_id]['gpr_by_label'] = GPR_D[GPR_id][
                            'gpr_by_label'
                        ].replace(gene_ids_sorted[x_], gene_labels_sorted[x_])
                        GPR_D[GPR_id]['gpr_by_name'] = GPR_D[GPR_id][
                            'gpr_by_name'
                        ].replace(gene_ids_sorted[x_], gene_names_sorted[x_])

                GPR_D[GPR_id]['gene_labels'] = gene_labels_sorted
                GPR_D[GPR_id]['miriam'] = None
                GPR_D[GPR_id]['annotation'] = {}
                GPR_D[GPR_id]['sbo'] = SBgpr.getSBOTermID()
                GPR_D[GPR_id]['notes'] = sbml_getNotes(SBgpr)
                if LOADANNOT:
                    if fbc_version < 3:
                        GPR_D[GPR_id]['annotation'] = sbml_readKeyValueDataAnnotation(SBgpr.getAnnotationString())
                    else:
                        # TODO bgoli deal with v3 extended annotation
                        GPRfbc = SBgpr.getPlugin("fbc")
                        GPR_D[GPR_id]['annotation'], GPR_D[GPR_id]['annotation_ext'] = sbml_readFBCv3KeyValuePairs(GPRfbc)

                    manot = sbml_getCVterms(SBgpr, model=False)

                    if manot != None:
                        GPR_D[GPR_id]['miriam'] = manot
                    del manot

        USE_NET_STOICH = False
        substrates = {}
        products = {}
        reagents = {}

        EXREAC = False
        reactionIDs.append(R_id)
        for sub in range(SBRe.getNumReactants()):
            spec = SBRe.getReactant(sub).getSpecies()
            stoi = -SBRe.getReactant(sub).getStoichiometry()
            if USE_NET_STOICH:
                if spec not in reagents:
                    reagents[spec] = float(stoi)
                else:
                    reagents[spec] += float(stoi)
            else:
                if spec not in substrates:
                    substrates[spec] = float(stoi)
                else:
                    substrates[spec] += float(stoi)

            if spec in boundary_species:
                EXREAC = True
        for pr in range(SBRe.getNumProducts()):
            spec2 = SBRe.getProduct(pr).getSpecies()
            stoi2 = SBRe.getProduct(pr).getStoichiometry()
            if USE_NET_STOICH:
                if spec2 not in reagents:
                    reagents[spec2] = float(stoi2)
                else:
                    reagents[spec2] += float(stoi2)
            else:
                if spec2 not in products:
                    products[spec2] = float(stoi2)
                else:
                    products[spec2] += float(stoi2)
            if spec2 in boundary_species:
                EXREAC = True
        R = CBModel.Reaction(SBRe.getId(), SBRe.getName(), reversible=SBRe.getReversible())
        reactionsReversability.append(SBRe.getReversible())

        if USE_NET_STOICH:
            for r in reagents:
                R.addReagent(
                    CBModel.Reagent('{}_{}'.format(SBRe.getId(), r), r, reagents[r])
                )
        else:
            if len(set(substrates).intersection(set(products))) == 0:
                substrates.update(products)
                products = {}
            for r in substrates:
                R.addReagent(
                    CBModel.Reagent('{}_{}'.format(SBRe.getId(), r), r, substrates[r])
                )
            for r in products:
                if r not in substrates:
                    R.addReagent(
                        CBModel.Reagent('{}_{}'.format(SBRe.getId(), r), r, products[r])
                    )
                else:
                    R.addReagent(
                        CBModel.Reagent(
                            '{}_{}_prod'.format(SBRe.getId(), r), r, products[r]
                        )
                    )
        del reagents, substrates, products

        if EXREAC:
            R.is_exchange = True
        R.annotation = {}
        if LOADANNOT:
            if fbc_version < 3:
                R.annotation = sbml_readKeyValueDataAnnotation(SBRe.getAnnotationString())
                # only dig for ancient annotation if not using V2
                if libsbml.getLibSBMLVersion() >= 51903:
                    node_txt = libsbml.XMLNode.convertXMLNodeToString(SBRe.getNotes())
                else:
                    node_txt = libsbml.XMLNode_convertXMLNodeToString(SBRe.getNotes())
                if fbc_version < 2 and R.annotation == {}:
                    R.annotation = sbml_readCOBRANote(node_txt)
            else:
                # TODO bgoli deal with v3 extended annotation
                #GPRfbc = SBgpr.getPlugin("fbc")
                R.annotation, R.annotation_ext = sbml_readFBCv3KeyValuePairs(RFBCplg)

            manot = sbml_getCVterms(SBRe, model=False)

            if manot != None:
                R.miriam = manot
            del manot
        setCBSBOterm(SBRe.getSBOTermID(), R)
        R.setNotes(sbml_getNotes(SBRe))
        if SBRe.getNumModifiers() > 0:
            for mo_ in range(SBRe.getNumModifiers()):
                R._modifiers_.append(SBRe.getModifier(mo_).getSpecies())
        REAC.append(R)

    if DEBUG:
        print('Reactions load: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # define compartments
    COMP = []
    for c_ in range(M.getNumCompartments()):
        SBcmp = M.getCompartment(c_)
        cid = SBcmp.getId()
        name = SBcmp.getName()
        size = SBcmp.getSize()
        if numpy.isnan(size) or size == None or size == '':
            size = None
        #volume = SBcmp.getVolume()
        dimensions = SBcmp.getSpatialDimensions()
        if dimensions == 0:
            print('Zero dimension compartment detected: {}'.format(cid))
            # zero dimension compartments make no sense and are assumed to be L2 artifacts
            dimensions = 3

        C = CBModel.Compartment(cid, name=name, size=size, dimensions=dimensions)

        if LOADANNOT:
            # TODO fix annotation
            if fbc_version < 3:
                C.annotation = sbml_readKeyValueDataAnnotation(SBcmp.getAnnotationString())
            else:
                # TODO bgoli deal with v3 extended annotation
                Cfbc = SBcmp.getPlugin("fbc")
                C.annotation, C.annotation_ext = sbml_readFBCv3KeyValuePairs(Cfbc)

            manot = sbml_getCVterms(SBcmp, model=False)

            if manot != None:
                C.miriam = manot
            del manot
        setCBSBOterm(SBcmp.getSBOTermID(), C)
        C.setNotes(sbml_getNotes(SBcmp))
        COMP.append(C)
        del cid, name, size, dimensions, C

    if DEBUG:
        print('Compartment load: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # extract fluxbounds
    if HAVE_FBC and fbc_version == 1:
        FB_data = []
        for fb_ in range(FBCplg.getNumFluxBounds()):
            SBFb = FBCplg.getFluxBound(fb_)
            fbd = {
                'reaction': SBFb.getReaction(),
                'operation': SBFb.getOperation(),
                'value': SBFb.getValue(),
            }
            fb_id = SBFb.getId()
            if fb_id not in ['', None]:
                fbd['id'] = fb_id
            FB_data.append(fbd)

        if DEBUG:
            print('FluxBounds load: {}'.format(round(time.time() - time0, 3)))
        time0 = time.time()

    # determine bound type and set default bounds for unbound reactions
    CONSTR = []
    boundReactionIDs = []
    LboundReactionIDs = []
    UboundReactionIDs = []
    AboundReactionIDs = []
    DefinedReactionIDs = []
    if HAVE_FBC and fbc_version == 1:
        cntr = 0
        for c in FB_data:
            if 'id' in c:
                newId = c['id']
                DefinedReactionIDs.append(c['reaction'])
            else:
                O = '{}'.format(cntr)
                if c['operation'] in ['greater', 'greaterEqual', '>', '>=']:
                    O = 'lower'
                    LboundReactionIDs.append(c['reaction'])
                elif c['operation'] in ['less', 'lessEqual', '<', '<=']:
                    O = 'upper'
                    UboundReactionIDs.append(c['reaction'])
                elif c['operation'] in ['equal', '=']:
                    O = 'equal'
                    AboundReactionIDs.append(c['reaction'])
                newId = '{}_{}_bnd'.format(c['reaction'], O)

            FB = CBModel.FluxBound(
                newId, c['reaction'], c['operation'], float(c['value'])
            )
            FB.setName(newId)
            CONSTR.append(FB)
            cntr += 1
            if c['reaction'] not in boundReactionIDs:
                boundReactionIDs.append(c['reaction'])
        ubcntr = 0
        for J in range(len(reactionIDs)):
            LBt = False
            UBt = False
            ABt = False
            if reactionIDs[J] in DefinedReactionIDs:
                LBt = UBt = ABt = True
            if reactionIDs[J] in LboundReactionIDs:
                LBt = True
            if reactionIDs[J] in UboundReactionIDs:
                UBt = True
            if reactionIDs[J] in AboundReactionIDs:
                ABt = True
            if not LBt and not UBt and not ABt:
                # print 'not LBt and not UBt and not ABt'
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'lower')
                if reactionsReversability[J]:
                    ##  print 'Adding reversible'
                    CONSTR.append(
                        CBModel.FluxBound(
                            newId, reactionIDs[J], 'greaterEqual', -numpy.inf
                        )
                    )
                else:
                    ##  print 'Adding irreversible'
                    CONSTR.append(
                        CBModel.FluxBound(newId, reactionIDs[J], 'greaterEqual', 0.0)
                    )
                ubcntr += 1
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'upper')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'lessEqual', numpy.inf)
                )
                ubcntr += 1
            elif not LBt and not ABt:
                # print 'not LBt and not ABt'
                # print reactionIDs[J]
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'lower')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'greaterEqual', -numpy.inf)
                )
                ubcntr += 1
                print('Added new lower bound', newId)
            elif not UBt and not ABt:
                # print 'not UBt and not ABt'
                # print reactionIDs[J]
                # newId = 'UC_%i' % ubcntr
                newId = '%s_%s_bnd' % (reactionIDs[J], 'upper')
                CONSTR.append(
                    CBModel.FluxBound(newId, reactionIDs[J], 'lessEqual', numpy.inf)
                )
                ubcntr += 1
                # print 'Added new upper bound', newId

    # TODO bgoli this has to be deprecated for FBC v 2+ by 0.9
    elif HAVE_FBC and fbc_version >= 2:
        timeFBV2 = time.time()
        for bnd in FB_data:
            FB = CBModel.FluxBound(bnd['id'], bnd['reaction'], bnd['operation'], bnd['value'])
            FB.annotation = bnd['annotation']
            FB.annotation_ext = bnd['annotation_ext']
            FB.miriam = bnd['miriam']
            FB.__param__ = bnd['parameter']
            FB.__sbo_term__ = bnd['sbo']
            FB.name = bnd['name']
            if FB.name == None or FB.name == '':
                FB.name = bnd['parameter']
            CONSTR.append(FB)
        print('FluxBounds process1: {}'.format(round(time.time() - timeFBV2, 3)))

    if DEBUG:
        print('FluxBounds process: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # Create parameters
    PARAM = []
    for p_ in PARAM_D:
        P = CBModel.Parameter(
            p_, PARAM_D[p_]['value'], PARAM_D[p_]['name'], PARAM_D[p_]['constant']
        )
        P.annotation = PARAM_D[p_]['annotation']
        P.annotation_ext = PARAM_D[p_]['annotation_ext']
        P.miriam = PARAM_D[p_]['miriam']
        P.__sbo_term__ = PARAM_D[p_]['sbo']
        P._associations_ = PARAM_D[p_]['association']
        P.setNotes(PARAM_D[p_]['notes'])
        if PARAM_D[p_]['is_fluxbound']:
            P._is_fluxbound_ = True
        PARAM.append(P)

    if DEBUG:
        print('Parameter process: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # GPRASSOC = {}
    if HAVE_FBC and LOADGENES:
        if fbc_version == 1 and __HAVE_FBA_ANOT_GENEASS__:
            SBGPR = FBCplg.getListOfGeneAssociations()
            for g_ in SBGPR:
                gprid = g_.getId()
                rid = g_.getReaction()
                GPR_D[gprid] = {}
                if g_.getAssociation() != None:
                    assoc = g_.getAssociation().toInfix()
                    if assoc == '' or assoc is None:
                        pass
                    else:
                        GPR_D[gprid]['reaction'] = rid
                        GPR_D[gprid]['gpr_by_id'] = assoc
                        GPR_D[gprid]['gpr_tree'] = sbml_getGPRasDictFBCv1(
                            ast.parse(assoc).body[0], {}
                        )

    # BUILD MODEL, we now need to do it here to link in the fluxobjectives defined in the objective functions
    if CUSTOM_MODEL_EXTENSION or CUSTOM_MODEL_METACLASS:
        fm = CBModelExtended(model_id)
    else:
        fm = CBModel.Model(model_id)
    if M.isSetMetaId():
        fm.__metaid__ = M.getMetaId()
    else:
        fm.setMetaId('meta_{}'.format(model_id))
    fm.name = model_name
    fm.setNotes(model_description)
    fm.annotation = sbml_readKeyValueDataAnnotation(M.getAnnotationString())
    fm.__FBC_STRICT__ = FBCstrict
    fm.__FBC_VERSION__ = fbc_version

    # try extract objective functions
    OBJFUNCout = []
    if HAVE_FBC:
        try:
            ACTIVE_OBJ = FBCplg.getActiveObjective().getId()
            print('INFO: Active objective:', ACTIVE_OBJ)
        except Exception as why:
            # AttributeError
            print('\nINFO: No active objective defined')
            print('\t', type(why))

        for of_ in range(FBCplg.getNumObjectives()):
            SBOf = FBCplg.getObjective(of_)
            OF = CBModel.Objective(SBOf.getId(), SBOf.getType())
            OF.setName(SBOf.getName())
            if OF.getId() == ACTIVE_OBJ:
                fm.addObjective(OF, active=True)
            else:
                fm.addObjective(OF, active=False)

            for ofl_ in range(SBOf.getNumFluxObjectives()):
                SBOfl = SBOf.getFluxObjective(ofl_)
                if SBOfl.getId() in [None, '']:
                    oid = '%s_%s_flobj' % (SBOf.getId(), SBOfl.getReaction())
                else:
                    oid = SBOf.getId()
                if fbc_version < 3:
                    Oflx = CBModel.FluxObjective(oid, SBOfl.getReaction(), float(SBOfl.getCoefficient()))
                else:
                    Oflx = CBModel.FluxObjective(oid, SBOfl.getReaction(), float(SBOfl.getCoefficient()), \
                                                 FBC3_VARIABLE_TYPES[SBOfl.getVariableType()])
                    print('vtype', Oflx.getType())
                Oflx.setName(SBOfl.getName())
                OF.addFluxObjective(Oflx)
            OBJFUNCout.append(OF)

        if DEBUG:
            print('ObjectiveFunction load: {}'.format(round(time.time() - time0, 3)))
        time0 = time.time()

    if DEBUG:
        print('GPR load: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    manot = sbml_getCVterms(M, model=True)
    if manot != None:
        fm.miriam = manot
    del manot

    fm.sourcefile = fname
    sbmh = M.getModelHistory()
    if sbmh != None:
        cdate = sbmh.getCreatedDate()
        if cdate != None:
            cdate = (
                cdate.getYear(),
                cdate.getMonth(),
                cdate.getDay(),
                cdate.getHour(),
                cdate.getMinute(),
                cdate.getSecond(),
            )
            fm.setCreatedDate(cdate)
        mdate = sbmh.getModifiedDate()
        if mdate != None:
            mdate = (
                mdate.getYear(),
                mdate.getMonth(),
                mdate.getDay(),
                mdate.getHour(),
                mdate.getMinute(),
                mdate.getSecond(),
            )
            fm.setModifiedDate(mdate)
        mCr = []
        for m_ in range(sbmh.getNumCreators()):
            sbc = sbmh.getCreator(m_)
            fm.addModelCreator(
                sbc.getGivenName(),
                sbc.getFamilyName(),
                sbc.getOrganisation(),
                sbc.getEmail(),
            )

    if DEBUG:
        print('Model build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    for co_ in COMP:
        fm.addCompartment(co_)
    if DEBUG:
        print('Compartment build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    for s_ in SPEC:
        fm.addSpecies(s_)
    if DEBUG:
        print('Species build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    for r_ in REAC:
        fm.addReaction(r_, create_default_bounds=False)
    if DEBUG:
        print('Reaction build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    fbexists = []
    timeFBV2 = time.time()
    for c_ in CONSTR:
        fm.addFluxBound(c_, fbexists=fbexists)
        fbexists.append("{}_{}".format(c_.getReactionId(), c_.getType()))
    print('FluxBounds process2: {}'.format(round(time.time() - timeFBV2, 3)))
    del fbexists
    if DEBUG:
        print('FluxBound build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    if DEBUG:
        print('Objective build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    for p_ in PARAM:
        try:
            fm.addParameter(p_)
        except RuntimeError:
            print('INFO: duplicate parameter id detected: {}'.format(p_.getId()))

    if DEBUG:
        print('Parameter build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()
    gene_labels = {}
    for g_ in GENE_D:
        gene_labels[GENE_D[g_]['id']] = GENE_D[g_]['label']

    if fbc_version == 1 and LOADGENES:
        for g_ in GPR_D:
            if 'reaction' in GPR_D[g_] and GPR_D[g_]['gpr_tree'] is not None:

                # fm.createGeneProteinAssociation(GPR_D[g_]['reaction'], GPR_D[g_]['gpr_by_id'], gid=g_,\
                # update_idx=False, altlabels=gene_labels)

                fm.createGeneProteinAssociationFromTree(
                    GPR_D[g_]['reaction'],
                    GPR_D[g_]['gpr_tree'],
                    gid=g_,
                    update_idx=False,
                    altlabels=gene_labels,
                )
                gpr = fm.getGPRassociation(g_)
                if gpr is not None:
                    gpr.setTree(GPR_D[g_]['gpr_tree'])
        fm.__updateGeneIdx__()
    elif fbc_version >= 2 and LOADGENES:
        # note we may want to add branches here for using indexes etc etc
        non_gpr_genes = []
        for g_ in GPR_D:
            if (
                'gpr_tree' in GPR_D[g_]
                and GPR_D[g_]['gpr_tree'] is not None
                and len(GPR_D[g_]['gpr_tree']) > 0
            ):
                # print(GPR_D[g_]['gpr_tree'])
                # fm.createGeneProteinAssociation(GPR_D[g_]['reaction'], GPR_D[g_]['gpr_by_id'], gid=g_, update_idx=False, altlabels=gene_labels)

                fm.createGeneProteinAssociationFromTree(
                    GPR_D[g_]['reaction'],
                    GPR_D[g_]['gpr_tree'],
                    gid=g_,
                    update_idx=False,
                    altlabels=gene_labels,
                )

                gpr = fm.getGPRassociation(g_)
                if gpr is not None:
                    gpr.annotation = GPR_D[g_]['annotation']
                    gpr.miriam = GPR_D[g_]['miriam']
                    gpr.__sbo_term__ = GPR_D[g_]['sbo']
                    gpr.setTree(GPR_D[g_]['gpr_tree'])
        fm.__updateGeneIdx__()
        for g_ in GENE_D:
            G = fm.getGene(g_)
            if G is not None:
                if G.getLabel() != GENE_D[g_]['label']:
                    print(G.getLabel(), GENE_D[g_]['label'])
                    G.setLabel(GENE_D[g_]['label'])
                G.name = GENE_D[g_]['name']
                G.annotation = GENE_D[g_]['annotation']
                G.__sbo_term__ = GENE_D[g_]['sbo']
                G.miriam = GENE_D[g_]['miriam']
                G.setNotes(GENE_D[g_]['notes'])
            else:
                if g_ is not None and g_ != 'None':
                    print(
                        'Gene {} is not part of a GPR association. Will create anyway!'.format(
                            g_
                        )
                    )
                    non_gpr_genes.append(g_)

        for ng_ in non_gpr_genes:
            G = CBModel.Gene(ng_, label=GENE_D[ng_]['label'], active=False)
            G.annotation = GENE_D[ng_]['annotation']
            G.__sbo_term__ = GENE_D[ng_]['sbo']
            G.miriam = GENE_D[ng_]['miriam']
            G.setNotes(GENE_D[ng_]['notes'])
            print(
                'WARNING: Non-gpr gene detected.', G.getId(), G.getLabel(), fm.getId()
            )

    if DEBUG:
        print('GPR build: {}'.format(round(time.time() - time0, 3)))
    time0 = time.time()

    # Groups support
    HAVE_GROUPS = False
    try:
        GRPplg = M.getPlugin(str("groups"))
        if GRPplg != None:
            HAVE_GROUPS = True
    except Exception as why:
        print(type(why))
        GRPplg = None
        HAVE_GROUPS = False

    if HAVE_GROUPS and GRPplg.getNumGroups() > 0:
        print('Groups support: {}'.format(GRPplg))
        print('Group.getNumGroups: {}'.format(GRPplg.getNumGroups()))

        s_ids = fm.getSpeciesIds()
        r_ids = fm.getReactionIds()
        c_ids = fm.getCompartmentIds()
        f_ids = fm.getFluxBoundIds()
        o_ids = fm.getObjectiveIds()
        g_ids = fm.getGeneIds()
        grp_ids = [GRPplg.getGroup(g).getId() for g in range(GRPplg.getNumGroups())]

        # print(grp_ids)
        sub_group_assoc = {}

        for g in range(GRPplg.getNumGroups()):
            GR = GRPplg.getGroup(g)
            grp_id = GR.getId()
            grp = CBModel.Group(grp_id)
            grp.setName(GR.getName())
            grp.setKind(GROUP_KINDS[GR.getKind()])
            if GR.getSBOTerm() != -1:
                grp.setSBOterm('SBO:{}'.format(str(GR.getSBOTerm()).zfill(7)))

            notes = GR.getNotesString()
            if notes != '' and notes != None:
                notes = notes.replace('<html:body>', '').replace('</html:body>', '')
                notes = notes.replace('<notes>', '').replace('</notes>', '').strip()
                grp.setNotes(notes)

            annostr = GR.getAnnotationString()
            if annostr != '' and annostr != None:
                grp.annotation = sbml_readKeyValueDataAnnotation(annostr)
            if DEBUG:
                print('\nid:', grp_id)
                print('name:', GR.getName())
                print('kind:', GROUP_KINDS[GR.getKind()])
                print('sboterm:', 'SBO:{}'.format(str(GR.getSBOTerm()).zfill(7)))
                print(
                    'annotations:',
                    sbml_readKeyValueDataAnnotation(GR.getAnnotationString()),
                )
                print('notes:', notes)

            LOM = GR.getListOfMembers()
            if LOM.getSBOTerm() != -1:
                grp.setSharedSBOterm('SBO:{}'.format(str(LOM.getSBOTerm()).zfill(7)))

            lom_notes = LOM.getNotesString()
            if lom_notes != '' and lom_notes != None:
                lom_notes = lom_notes.replace('<html:body>', '').replace(
                    '</html:body>', ''
                )
                lom_notes = (
                    lom_notes.replace('<notes>', '').replace('</notes>', '').strip()
                )
                grp.setSharedNotes(lom_notes)

            annostr = LOM.getAnnotationString()
            if annostr != '' and annostr != None:
                grp._member_attributes_.annotation = sbml_readKeyValueDataAnnotation(
                    annostr
                )

            if DEBUG:
                print(
                    'sboterm (shared):', 'SBO:{}'.format(str(LOM.getSBOTerm()).zfill(7))
                )
                print(
                    'annotations (shared):',
                    sbml_readKeyValueDataAnnotation(LOM.getAnnotationString()),
                )
                print('notes (shared):', lom_notes)
                print('num members:', GR.getNumMembers())

            for m in range(GR.getNumMembers()):
                idr = GR.getMember(m).getIdRef()
                if DEBUG:
                    print('\tmember idRef:', GR.getMember(m).getIdRef())
                if idr in s_ids:
                    grp.addMember(fm.species[s_ids.index(idr)])
                elif idr in r_ids:
                    grp.addMember(fm.reactions[r_ids.index(idr)])
                elif idr in f_ids:
                    grp.addMember(fm.flux_bounds[f_ids.index(idr)])
                elif idr in g_ids:
                    grp.addMember(fm.genes[g_ids.index(idr)])
                elif idr in c_ids:
                    grp.addMember(fm.compartments[c_ids.index(idr)])
                elif idr in o_ids:
                    grp.addMember(fm.objectives[o_ids.index(idr)])
                elif idr in grp_ids:
                    # print(idr)
                    # print(grp_id)
                    # scanning for subgroup members to assign later
                    if grp_id not in sub_group_assoc:
                        sub_group_assoc[grp_id] = [idr]
                    else:
                        sub_group_assoc[grp_id].append(idr)
                else:
                    print(
                        'Skipping group \"{}\" member \"{}\", it is an incompatible type.'.format(
                            GR.getId(), idr
                        )
                    )

            fm.addGroup(grp)
        del s_ids, r_ids, c_ids, f_ids, o_ids, g_ids

        # adding subgroups
        for grp_ in sub_group_assoc:
            GRP = fm.getGroup(grp_)
            for sg_ in sub_group_assoc[grp_]:
                GRP.addMember(fm.getGroup(sg_))
                # print('Adding: {} to {}'.format(sg_, grp_))
        # print(sub_group_assoc)
        del sub_group_assoc

    print('')
    if len(CONSTR) < 1:
        print('INFO: No FBC flux bounds were defined.')
        # time.sleep(1)
    if len(OBJFUNCout) < 1:
        print('INFO: No FBC objective functions were defined.')
        # time.sleep(1)
    fm._SBML_LEVEL_ = 3

    del SPEC, GENE_D, REAC, PARAM, PARAM_D
    try:
        if NMATRIX_TYPE is not None:
            fm.buildStoichMatrix(matrix_type=NMATRIX_TYPE)
    except Exception as why:
        print(type(why))
        print('\nINFO: unable to construct stoichiometric matrix')
    if DEBUG:
        print('Nmatrix build: {}'.format(round(time.time() - time0, 3)))

    print('SBML3 load time: {}\n'.format(round(time.time() - time00, 3)))

    # del M, D

    if not return_sbml_model:
        return fm
    else:
        print('SBML object return disabled')
        return fm, None


def sbml_readCOBRANote(s):
    """
    Parses a COBRA style note from a XML string

     - *s* an XML string

    """

    new_ann = {}
    cobra_p = re.findall(re_html_cobra_p, s)
    html_c_p = re.findall(re_html_p, s)
    if len(cobra_p) > 0:
        if __DEBUG__:
            print(cobra_p)
        for p in cobra_p:
            ps = (
                p.replace('<html:p>', '')
                .replace('</html:p>', '')
                .replace('&lt;', '<')
                .replace('&gt;', '>')
                .split(':', 1)
            )
            if len(ps) == 2:
                new_ann.update({ps[0].strip(): ps[1].strip()})
    elif len(html_c_p) > 0:
        if __DEBUG__:
            print(cobra_p)
        for p in html_c_p:
            ps = (
                p.replace('<p>', '')
                .replace('</p>', '')
                .replace('&lt;', '<')
                .replace('&gt;', '>')
                .split(':', 1)
            )
            if len(ps) == 2:
                new_ann.update({ps[0].strip(): ps[1].strip()})
    elif '<span xmlns="http://www.w3.org/1999/xhtml">' in s:
        hPs = re.findall(re_html_p, s.replace('\n', ''))
        if __DEBUG__:
            print(hPs)
        for p in hPs:
            ps = re.findall(re_html_span, p)
            ps = [
                p.replace('<span>', '')
                .replace('</span>', '')
                .replace('&lt;', '<')
                .replace('&gt;', '>')
                .strip()
                for p in ps
            ]
            if len(ps) == 2:
                new_ann.update({ps[0]: ps[1]})
    # print new_ann
    return new_ann


def sbml_getCVterms(sb, model=False):
    """
    Get the MIRIAM compliant CV terms and return a MIRIAMAnnotation or None

     - *sb* a libSBML SBase derived object
     - *model* is this a BQmodel term

    """
    out = {}
    if sb.getNumCVTerms() > 0:
        for c_ in range(sb.getNumCVTerms()):
            cvt = sb.getCVTerm(c_)
            resrc = []
            if cvt.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER:
                qual = BQ2CBMMAP[cvt.getBiologicalQualifierType()]
            else:
                qual = BQM2CBMMAP[cvt.getModelQualifierType()]
            for r_ in range(cvt.getNumResources()):
                uri = cvt.getResourceURI(r_)
                if qual != None and uri != None:
                    resrc.append(uri)
            out[qual] = tuple(resrc)
    if len(out) > 0:
        man = MIRIAMannotation()
        for q_ in out:
            for r_ in out[q_]:
                man.addIDorgURI(q_, r_)
        return man
    else:
        return None


def sbml_setCVterms(sb, uridict, model=False):
    """
    Add MIRIAM compliant CV terms to a sbml object from a CBM object

     - *sb* a libSBML SBase derived object
     - *uridict* a dictionary of uri's as produced by getAllMIRIAMUris()
     - *model* is this a BQmodel term [deprecated attribute, ignored and autodetected]

    """
    for t_ in uridict:
        if len(uridict[t_]) > 0:
            global cv
            # print(t_)
            # print(uridict[t_])
            if t_ in CBM2BQMAP:
                cv = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
                cv.setBiologicalQualifierType(CBM2BQMAP[str(t_)])
                # print(CBM2BQMAP[t_])
            else:
                cv = libsbml.CVTerm(libsbml.MODEL_QUALIFIER)
                cv.setModelQualifierType(CBM2BQMMAP[str(t_)])
            for u_ in uridict[t_]:
                # print(u_)
                # print(uridict[t_])
                if cv.addResource(str(u_)) != libsbml.LIBSBML_OPERATION_SUCCESS:
                    print(
                        'INFO failure adding MIRIAM resource {} to: {}'.format(
                            u_, sb.getId()
                        )
                    )

            if sb.addCVTerm(cv) != libsbml.LIBSBML_OPERATION_SUCCESS:
                print('INFO: failure adding MIRIAM term: {}'.format(sb.addCVTerm(cv)))
