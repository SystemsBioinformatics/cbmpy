"""
CBMPy: CBModel module
=====================
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
Last edit: $Author: bgoli $ ($Id: CBModel.py 696 2019-07-29 21:59:43Z bgoli $)

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
#from __future__ import unicode_literals

import numpy, re, time, weakref, copy, json, ast, os

try:
    from urllib2 import quote as urlquote
    from urllib2 import unquote as urlunquote
except ImportError:
    from urllib.parse import quote as urlquote
    from urllib.parse import unquote as urlunquote

global GENE_CNTR
GENE_CNTR = 0

try:
    import cPickle as pickle
except ImportError:
    import pickle


HAVE_SYMPY = False
try:
    import sympy
    if int(sympy.__version__.split('.')[1]) >= 7 and int(sympy.__version__.split('.')[2]) >= 5:
        HAVE_SYMPY = True
    elif int(sympy.__version__.split('.')[0]) >= 1:
        HAVE_SYMPY = True
    else:
        del sympy
        print('\nWARNING: SymPy version 0.7.5 or newer is required for symbolic matrix support.')
except ImportError:
    HAVE_SYMPY = False
    print('SymPy not install (only required for optional, symbolic matrix support).')

HAVE_SCIPY = False
try:
    import scipy
    from scipy.sparse import csr_matrix
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


#from .CBDataStruct import (StructMatrixLP, MIRIAMannotation, MIRIAMModelAnnotation)
from .CBDataStruct import (StructMatrixLP, MIRIAMannotation)
from .CBCommon import (checkChemFormula, extractGeneIdsFromString, getGPRasDictFromString,\
                       binHash, fixId, checkId, createAssociationDictFromNode)

from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']


class Fbase(object):
    """
    Base class for CB Model objects

    """
    id = None
    name = None
    annotation = None
    compartment = None
    miriam = None
    notes = ''
    __TRASH__ = None
    __objref__ = None
    __metaid__ = None
    __sbo_term__ = None
    __text_encoding__ = 'utf8'

    __ENABLE_GLOBAL_WEAKREF__ = True
    ##  __pre__ = ''
    ##  __post__ = ''


    def __getstate__(self):
        """
        Internal method that should allow our weakrefs to be 'picklable'

        # overloaded by Model

        """

        if '__objref__' not in self.__dict__:
            return self.__dict__
        else:
            cpy = self.__dict__.copy()
            cpy['__objref__'] = None
            return cpy

    def __setObjRef__(self, o):
        """
        Sets a weakref to another object

         - *o* the referred to object

        """
        self.__objref__ = weakref.ref(o)

    def __unsetObjRef__(self):
        """
        UnSets a weakref to another object

         - *o* the referred to object

        """
        self.__objref__ = None

    def getPid(self):
        """
        Return the object ID.

        """
        return self.getId()

    def getId(self):
        """
        Return the object ID.

        """
        return self.id

    def getMetaId(self):
        """
        Return the object metaId.

        """
        if self.__metaid__ == None:
            self.__metaid__ = 'meta_{}'.format(self.id)
        return self.__metaid__

    def getName(self):
        """
        Return the object name.

        """
        return self.name

    def getNotes(self):
        """
        Return the object's notes

        """
        #return self.__urlDecode(self.notes)
        return self.notes

    def getAnnotations(self):
        """
        Return the object annotation dictionary

        """
        return self.annotation

    def getAnnotation(self, key):
        """
        Return the object annotation associated with:

         - *key* the annotation key

        """
        if key in self.annotation:
            return self.annotation[key]
        else:
            return None

    def hasAnnotation(self, key):
        """
        Returns a boolean representing the presence/absence of the key in the objext annotation

         - *key* the annotation key

        """
        if key in self.annotation:
            return True
        else:
            return False

    def setName(self, name):
        """
        Set the object name:

         - *name* the name string

        """
        self.name = str(name)

    def setNotes(self, notes):
        """
        Sets the object's notes:

         - *notes* the note string, should preferably be (X)HTML for SBML

        """
        #self.notes = self.__urlEncode(notes)
        #try:
            ##self.notes = notes.decode(errors='ignore')
            #self.notes = notes
        #except AttributeError:
            #self.notes=notes
        self.notes = notes

    def setAnnotation(self, key, value):
        """
        Set an objects annotation as a key : value pair.

         - *key* the annotation key
         - *value* the annotation value

        """
        assert self.annotation != None, '\nThis class has no annotation field'
        self.annotation.update({key : value})
        #if type(value) != list:
            #self.annotation.update({key : value})
        #elif type(value) == list:
            #if len(value) == 1:
                #self.annotation.update({key : value[0]})
            #else:
                #value = ['&apos;{}&apos;'.format(str(v).strip()) for v in value]
                #self.annotation.update({key : '['+', '.join(map(str, value))+']'})
                ##self.annotation.update({key : str(json.dumps(value))})
                ##self.annotation.update({key : value})

    def deleteAnnotation(self, key):
        """
        Unsets (deltes) an objects annotation with key

         - *key* the annotation key

        """
        assert key != self.annotation, '\nAnnotation key {} does not exist'.format(key)
        self.annotation.pop(key)

    def __checkId__(self, cid=None):
        """
        Checks the validity of the object id unless cid is provided, in which case it checks the provided id

         - *cid* [default=None] an optional Id to test for validity, if None then the object id is used

        """
        if cid is None:
            cid = self.id

        try:
            exec('{} = 0'.format(cid))
        except SyntaxError:
            #print('ERROR: Syntax error, \"{}\" is an invalid object identifier'.format(cid))
            return False

        ## this was the old way of doing it
        #cntr = 0
        #for c in cid:
            #if cntr == 0 and c.isalpha() or c == '_':
                #pass
            #elif cntr > 0 and c.isalnum() or c == '_':
                #pass
            #else:
                #print('\"{}\" is an invalid character in id: \"{}\"'.format(c, cid))
                ##raise RuntimeWarning('\n\"{}\" is an invalid character in id: \"{}\"'.format(c, cid))
                #return False
            #cntr += 1
        return True

    def setPid(self, fid):
        """
        Sets the object Id

         - *fid* a valid c variable style id string

        """
        self.setId(fid)

    def setId(self, fid):
        """
        Sets the object Id

         - *fid* a valid c variable style id string


         Reimplemented by @Reaction, @Species, @Compartment, @Gene

        """
        fid = str(fid)
        if fid == self.id:
            return

        if not self.__checkId__(fid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(fid))

        if self.__objref__ is not None:
            if fid not in self.__objref__().__global_id__:
                self.__objref__().__changeGlobalId__(self.id, fid, self)
                self.id = fid
            else:
                print('ERROR: setId() - object with id \"{}\" already exists ... ID *not* set.'.format(fid))
        else:
            self.id = fid

    def setMetaId(self, mid=None):
        """
        Sets the object Id

         - *mid* [default=None] a valid c variable style metaid string, if None it will be set as meta+id

        """
        if mid == None:
            self.__metaid__ = 'meta_{}'.format(self.id)
        else:
            self.__metaid__ = str(mid)

    def clone(self):
        """
        Return a clone of this object. Cloning performs a deepcop on the object which will also clone
        any objects that exist as attributes of this object, in other words an independent copy of the
        original. If this is not the desired behaviour override this method when subclassing or implement
        your own.

        """
        if type(self.__TRASH__) == dict:
            self.__TRASH__.clear()
        else:
            self.__TRASH__ = None
        return copy.deepcopy(self)

    def serialize(self, protocol=0):
        """
        Serialize object, returns a string by default

         - *protocol* [default=0] serialize to a string or binary if required,
                      see pickle module documentation for details

        # Reimplemented in Model

        """
        return pickle.dumps(self, protocol=protocol)

    def serializeToDisk(self, filename, protocol=2):
        """
        Serialize to disk using pickle protocol:

         - *filename* the name of the output file
         - *protocol* [default=2] serialize to a string or binary if required,
                      see pickle module documentation for details

        # Reimplemented in Model

        """
        F = open(filename, 'wb')
        pickle.dump(self, F, protocol=protocol)
        F.close()

    def addMIRIAMannotation(self, qual, entity, mid):
        """
        Add a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158 or fully qualifies url (if only_qual_uri)

        """
        if self.miriam == None:
            self.miriam = MIRIAMannotation()
        self.miriam.addMIRIAMannotation(qual, entity, mid)

    def addMIRIAMuri(self, qual, uri):
        """
        Add a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *uri* the fully qualified entity id e.g. http://identifiers.org/chebi/CHEBI:12345 (no validity checking is done)

        """
        if self.miriam == None:
            self.miriam = MIRIAMannotation()
        self.miriam.addIDorgURI(qual, uri)

    def getMIRIAMannotations(self):
        """
        Returns a dictionary of all MIRIAM annotations associated with this object
        or None of there are none defined.

        """
        if self.miriam != None:
            try:
                return self.miriam.getAllMIRIAMUris()
            except Exception as ex:
                print(ex)
                return None
        else:
            return None

    def deleteMIRIAMannotation(self, qual, entity, mid):
        """
        Deletes a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158

        """
        self.miriam.deleteMIRIAMannotation(qual, entity, mid)

    def getCompartmentId(self):
        """
        Return the compartment id where this element is located

        """
        return self.compartment

    def setCompartmentId(self, compartment):
        """
        Set the compartment id where this element is located

        """
        self.compartment = compartment

    def getSBOterm(self):
        """
        Return the SBO term for this object.

        """
        return self.__sbo_term__

    def setSBOterm(self, sbo):
        """
        Set the SBO term for this object.

         - *sbo* the SBOterm with format: SBO:nnnnnnn"

        """
        assert sbo.startswith('SBO:') and len(sbo.split(':')[1]) == 7, 'SBOterm must have the form: SBO:nnnnnnn'
        self.__sbo_term__ = sbo


    def __urlEncode(self, txt):
        """
        URL encodes a string.

        """
        try:
            txt = urlquote(txt.encode(self.__text_encoding__, errors='replace'), safe='')
        except UnicodeDecodeError as why:
            pass
            #print(txt)
        return txt

    def __urlDecode(self, txt):
        """
        Decodes a URL encoded string

        """
        return urlunquote(txt)

class Model(Fbase):
    """
    Container for constraint based model, adds methods for manipulating:

     - objectives
     - constraints
     - reactions
     - species
     - compartments
     - groups
     - parameters
     - N a structmatrix object

    """

    prefix = ''
    suffix = ''
    objectives = None
    obj_func = None
    activeObjIdx = None
    flux_bounds = None
    reactions = None
    species = None
    compartments = None
    groups = None
    genes = None
    __genes_idx__ = None
    __single_gene_effect_map__ = None
    gpr = None
    parameters = None
    N = None
    sourcefile = ''
    description = ''
    user_constraints = None
    CM = None
    sensitivity = None
    ##  optValue = None
    ##  optType = None
    SCALED_REDUCED_COSTS = False
    SOLUTION_STATUS_INT = None
    SOLUTION_STATUS = None
    DATE_CREATED = None
    DATE_MODIFIED = None
    MODEL_CREATORS = None
    __gene_deactivated_reactions__ = None
    _SBML_LEVEL_ = None
    __FBC_VERSION__ = 1
    __FBC_STRICT__ = True
    #__objref__ = None
    __global_id__ = None
    __modified__ = False
    __check_gene_activity__ = False
    __CUSTOM_MODEL_EXTENSION__ = None
    __CUSTOM_MODEL_METACLASS__ = None


    def __init__(self, pid):
        """
        Container for constraint based model.

         - *pid* a string with the model ID

        Adds methods for manipulating:

         - objectives
         - constraints
         - reactions
         - species
         - compartments
         - groups
         - N a structmatrix object

        """
        pid = str(pid)
        self.setId(pid)
        self.objectives = []
        self.flux_bounds = []
        self.reactions = []
        self.species = []
        self.compartments = []
        self.groups = []
        self.genes = []
        self.__genes_idx__ = []
        self.gpr = []
        self.parameters = []
        self.annotation = {}
        self.__TRASH__ = {}
        self.MODEL_CREATORS = {}
        self.__gene_deactivated_reactions__ = {}
        self.compartment = None
        self.__delattr__('compartment')
        self.__setGlobalIdStore__()

    def __setGlobalIdStore__(self):
        """
        Does exactly what the function name says (creates/replaces the globalId store)

        """
        # testing weakref dictionary, needs work ...
        if self.__ENABLE_GLOBAL_WEAKREF__:
            self.__global_id__ = weakref.WeakValueDictionary({self.getId() : self})
        else:
            self.__global_id__ = {self.getId() : None}

    def clone(self):
        """
        Return a clone of this object.

        """
        tzero = time.time()
        if type(self.__TRASH__) == dict:
            self.__TRASH__.clear()
        else:
            self.__TRASH__ = None

        cpy = copy.deepcopy(self)
        cpy.__global_id__.clear()
        cpy.__setModelSelf__()
        cpy.__setGlobalIdStore__()
        cpy.__populateGlobalIdStore__()
        # TODO: try make the global ID store work exclusively with __setState__ and __getState__ for both clone/serialise
        self.__setGlobalIdStore__()
        self.__populateGlobalIdStore__()

        print('Model clone time: {}'.format(time.time()-tzero))
        return cpy

    def __populateGlobalIdStore__(self):
        """
        This method populates the globalID store after a clone This is a utility function that comes into play when a model is cloned.

        NB: synch with __setModelSelf__()

        """
        for r in self.reactions:
            r.__setObjRef__(self)
            self.__global_id__[r.getId()] = r
            for rr in r.reagents:
                self.__global_id__[rr.getId()] = rr
        for s in self.species:
            self.__global_id__[s.getId()] = s
        for fb in self.flux_bounds:
            self.__global_id__[fb.getId()] = fb
        for o in self.objectives:
            self.__global_id__[o.getId()] = o
            for fo in o.fluxObjectives:
                self.__global_id__[fo.getId()] = fo
        for c in self.compartments:
            self.__global_id__[c.getId()] = c
        for gp in self.gpr:
            self.__global_id__[gp.getId()] = gp
        for g in self.genes:
            self.__global_id__[g.getId()] = g
        for p in self.parameters:
            self.__global_id__[p.getId()] = p

    def __setModelSelf__(self):
        """
        This method sets the model reference (updates the weakref) to the current instance. This is a
        utility function that mostly comes into play when a model is cloned or objects are mixed between models.

        NB: synch with __unsetModelSelf__() and __populateGlobalIdStore__()

        """
        for r in self.reactions:
            r.__setObjRef__(self)
            for rr in r.reagents:
                rr.__setObjRef__(self)
        for s in self.species:
            s.__setObjRef__(self)
        for fb in self.flux_bounds:
            fb.__setObjRef__(self)
        for o in self.objectives:
            o.__setObjRef__(self)
            for fo in o.fluxObjectives:
                fo.__setObjRef__(self)
        for c in self.compartments:
            c.__setObjRef__(self)
        for gp in self.gpr:
            gp.__setObjRef__(self)
        for g in self.genes:
            g.__setObjRef__(self)
        for p in self.parameters:
            p.__setObjRef__(self)


    def __unsetModelSelf__(self):
        """
        This method unsets the model reference (deletes the weakref). This is a
        utility function that mostly comes into play when a model is cloned or objects are mixed between models.

        NB: synch with __setModelSelf__() and __populateGlobalIdStore__()

        """
        for r in self.reactions:
            r.__unsetObjRef__()
            for rr in r.reagents:
                rr.__unsetObjRef__()
        for s in self.species:
            s.__unsetObjRef__()
        for fb in self.flux_bounds:
            fb.__unsetObjRef__()
        for o in self.objectives:
            o.__unsetObjRef__()
            for fo in o.fluxObjectives:
                fo.__unsetObjRef__()
        for c in self.compartments:
            c.__unsetObjRef__()
        for gp in self.gpr:
            gp.__unsetObjRef__()
        for g in self.genes:
            g.__unsetObjRef__()
        for p in self.parameters:
            p.__unsetObjRef__()

    def __getstate__(self):
        """
        Internal method that should allow our weakrefs to be 'picklable'

        # overloaded by Model

        """

        self.__global_id__ = None
        if '__objref__' not in self.__dict__:
            return self.__dict__
        else:
            cpy = self.__dict__.copy()
            cpy['__objref__'] = None
            return cpy

    def __setstate__(self, dic):
        """
        Internal method that allows our weakrefs to be 'picklable'

        """
        self.__dict__ = dic
        self.__setModelSelf__()
        self.__setGlobalIdStore__()
        self.__populateGlobalIdStore__()

    def serialize(self, protocol=0):
        """
        Serialize object, returns a string by default

         - *protocol* [default=0] serialize to a string or binary if required,
                      see pickle module documentation for details

        # overloaded in CBModel

        """
        s = pickle.dumps(self, protocol=protocol)
        self.__setGlobalIdStore__()
        self.__populateGlobalIdStore__()
        return s

    def serializeToDisk(self, filename, protocol=2):
        """
        Serialize to disk using pickle protocol:

         - *filename* the name of the output file
         - *protocol* [default=2] serialize to a string or binary if required,
                      see pickle module documentation for details

        # overloaded in CBModel

        """
        F = open(filename, 'wb')
        pickle.dump(self, F, protocol=protocol)
        F.close()
        self.__setGlobalIdStore__()
        self.__populateGlobalIdStore__()

    def addMIRIAMannotation(self, qual, entity, mid):
        """
        Add a qualified MIRIAM annotation or entity:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158

        """
        if self.miriam == None:
            self.miriam = MIRIAMannotation()
        self.miriam.addMIRIAMannotation(qual, entity, mid)

    def setDescription(self, html):
        """
        Sets the model description which translates into the SBML <notes> field.

         - *html* any valid html or the empty string to clear ''

        """
        self.description = html
        #self.notes = self.__urlEncode(html)

    def getDescription(self):
        """
        Returns the model description which was stored in the SBML <notes> field

        """

        #return self.__urlDecode(self.notes)
        return self.notes

    def setCreatedDate(self, date=None):
        """
        Set the model created date tuple(year, month, day, hour, minute, second)

        - *date* [default=None] default is now (automatic) otherwise (year, month, day, hour, minute, second) e.g. (2012, 09, 24, 13, 34, 00)

        """
        if date == None:
            lt = time.gmtime()
            self.DATE_CREATED = (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)
        else:
            assert len(date) == 6, '\nInvalid date'
            self.DATE_CREATED = tuple([int(abs(t)) for t in date])
        if self.DATE_MODIFIED == None:
            self.setModifiedDate()

    def setModifiedDate(self, date=None):
        """
        Set the model modification date: tuple(year, month, day, hour, minute, second)

        - *date* [default=None] default is now (automatic) otherwise (year, month, day, hour, minute, second) e.g. (2012, 09, 24, 13, 34, 00)

        """
        if date == None:
            lt = time.gmtime()
            self.DATE_MODIFIED = (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec)
        else:
            assert len(date) == 6, '\nInvalid date'
            self.DATE_MODIFIED = tuple([int(abs(t)) for t in date])
        if self.DATE_CREATED == None:
            self.setCreatedDate()


    def addModelCreator(self, firstname, lastname, organisation=None, email=None):
        """
        Add a model creator to the list of model creators, only the first and fmaily names are mandatory:

         - *firstname*
         - *lastname*
         - *organisation* [default=None]
         - *email*  [default=None]

        """
        self.MODEL_CREATORS.update({firstname.title()+lastname.title() : {
            'firstname' : firstname,
            'lastname' : lastname,
            'organisation' : organisation,
            'email' : email
        }
                                    })

    def getModelCreators(self):
        """
        Return model creator information

        """
        out = []
        mc = self.MODEL_CREATORS
        for c_ in mc:
            out.append((mc[c_]['firstname'], mc[c_]['lastname'], mc[c_]['organisation'], mc[c_]['email']))
        return tuple(out)

    def addObjective(self, obj, active=False):
        """
        Add an instantiated Objective object to the FBA model

        - *obj* an instance of the Objective class
        - *active* [default=False] flag this objective as the active objective (fba.activeObjIdx)

        """
        assert type(obj) == Objective, '\nERROR: requires an Objective object, not something of type {}'.format(type(obj))
        assert obj.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(obj.__objref__).split('to')[1][1:-1])
        print('Adding objective: {}'.format(obj.id))
        obj.__objref__ = weakref.ref(self)
        if obj.getId() in self.__global_id__:
            raise RuntimeError('Duplicate obj ID detected: {}'.format(obj.getId()))
        else:
            self.__pushGlobalId__(obj.getId(), obj)
        self.objectives.append(obj)
        self.obj_func = self.objectives[-1]
        if active:
            self.activeObjIdx = len(self.objectives)-1

    def createObjectiveFunction(self, rid, coefficient=1, osense='maximize', active=True, delete_current_obj=True):
        """
        Create a single variable objective function:

         - **rid** The
         - **coefficient** [default=1]
         - **osense** [default='maximize']
         - **active** [default=True]
         - **delete_current_obj** [default=True]

        """
        new_obj_id = rid + '_objective'
        assert new_obj_id not in self.getObjectiveIds(), '\nObjective Id %s already exists' % (new_obj_id)
        if delete_current_obj:
            for o in self.getObjectiveIds():
                self.deleteObjective(o)
        obj = Objective(new_obj_id, osense)
        FO = FluxObjective('{}_{}_fluxobj'.format(new_obj_id, rid), rid, coefficient)
        self.addObjective(obj, active=active)
        obj.addFluxObjective(FO)

    def createSpecies(self, sid, boundary=False, name='', value=float('nan'), compartment=None, charge=None, chemFormula=None):
        """
        Create a new species and add it to the model:

         - **id** the unique species id
         - **boundary** [default=False] whether the species is a variable (False) or is a boundary parameter (fixed)
         - **name** [default=''] the species name
         - **value** [default=nan] the value *not currently used*
         - **compartment** [default=None] the compartment the species is located in
         - **charge** [default=None] the species charge
         - **chemFormula** [default=None] the chemical formula

        """

        assert sid not in self.getSpeciesIds(), '\nSpecies id\'s must be unique: {} exists'.format(sid)

        S = Species(sid, boundary, name, value, compartment, charge, chemFormula)
        self.addSpecies(S)


    def createCompartment(self, cid, name=None, size=1, dimensions=3, volume=None):
        """
        Create a new compartment and add it to the model if the id does not exist

        - *cid* compartment id
        - *name* [None] compartment name
        - *size* [1] compartment size
        - *dimensions* [3] compartment size dimensions
        - *volume* [None] compartment volume

        """
        if cid not in self.getCompartmentIds():
            C = Compartment(cid, name, size, dimensions, volume)
            self.addCompartment(C)
        else:
            print('Error: compartment id \"{}\"'.format(cid))



    def createReaction(self, rid, name=None, reversible=True, create_default_bounds=True, silent=False):
        """
        Create a new blank reaction and add it to the model:

         - *id* the unique reaction ID
         - *name* the reaction name
         - *reversible* [default=True] the reaction reversibility. True is reversible, False is irreversible
         - *create_default_bounds* create default reaction bounds, irreversible 0 <= J <= INF, reversable -INF <= J <= INF
         - *silent* [default=False] if enabled this disables the printing of information messages

        """

        assert rid not in self.getReactionIds(), '\nReaction ID %s already exists' % rid
        self.addReaction(Reaction(rid, name, reversible), create_default_bounds=create_default_bounds, silent=silent)
        if not silent:
            print('Add reagents with cmod.createReactionReagent({}, metabolite, coefficient)'.format(rid))

    def createReactionReagent(self, reaction, metabolite, coefficient, silent=False):
        """
        Add a reagent to an existing reaction, both reaction and metabolites must exist

         - *reaction* a reaction id
         - *metabolite* a species/metabolite id
         - *coefficient* the reagent coefficient

        """
        R = self.getReaction(reaction)
        S = self.getSpecies(metabolite)
        assert R != None, '\nReaction {} does not exist'.format(reaction)
        assert S != None, '\nMetabolite {} does not exist'.format(metabolite)
        R.createReagent(metabolite, coefficient)
        # this is handled by dynamic code
        #S.setReagentOf(reaction)
        if not silent:
            print('Added \"{}\" as a reagent of reaction \"{}\"'.format(metabolite, reaction))

    def createReactionLowerBound(self, reaction, value):
        """
        Create a new lower bound for a reaction: value <= reaction

         - **reaction** the reaction id
         - **value** the value of the bound

        """

        bnds = self.getReactionBounds(reaction)
        assert bnds[1] == None and bnds[3] == None, '\nLower or equality bound exists for reaction: %s' % reaction
        newId = '%s_%s_bnd'% (reaction, 'lower')
        self.addFluxBound(FluxBound(newId, reaction, 'greaterEqual', value))

    def createReactionUpperBound(self, reaction, value):
        """
        Create a new upper bound for a reaction: reaction <= value

         - **reaction** the reaction id
         - **value** the value of the bound

        """

        bnds = self.getReactionBounds(reaction)
        assert bnds[2] == None and bnds[3] == None, '\nUpper or equality bound exists for reaction: %s' % reaction
        newId = '%s_%s_bnd'% (reaction, 'upper')
        self.addFluxBound(FluxBound(newId, reaction, 'lessEqual', value))

    def createReactionBounds(self, reaction, lb_value, ub_value):
        """
        Create a new lower bound for a reaction: value <= reaction

         - **reaction** the reaction id
         - **lb_value** the value of the lower bound
         - **ub_value** the value of the upper bound

        """
        bnds = self.getReactionBounds(reaction)
        assert bnds[1] == None, '\nLower bound exists for reaction: %s' % reaction
        assert bnds[2] == None, '\nUpper bound exists for reaction: %s' % reaction
        assert bnds[3] == None, '\nEquality bound exists for reaction: %s' % reaction
        assert lb_value <= ub_value, '\nLower bound must be less than or equal to upper bound.'
        newId = '%s_%s_bnd'% (reaction, 'lower')
        self.addFluxBound(FluxBound(newId, reaction, 'greaterEqual', lb_value))
        newId = '%s_%s_bnd'% (reaction, 'upper')
        self.addFluxBound(FluxBound(newId, reaction, 'lessEqual', ub_value))

    def addFluxBound(self, fluxbound, fbexists=None):
        """
        Add an instantiated FluxBound object to the FBA model

        - *fluxbound* an instance of the FluxBound class
        - *fbexists* [default=None] this is a list of strings which contains fluxbounds that have been added to the model, see sample code below.
        The format of the string is 'reactionid_boundtype'

        ```python
        fbexists = []
        for fluxbound in list_of_fluxbounds:
            model.addFluxBound(fluxbound, fbexists=fbexists)
            fbexists.append("{}_{}".format(fluxbound.getReactionId(), fluxbound.getType()))
        ```

        Using the fbexists list drastically reduces the time it takes to add fluxbounds but circumvents any sort of existence checking and should *only* be
        used as shown above when constructing a model from scratch. Alternativel, you need to prepopulate fbexists with existing fluxbound component information:

        ```python
        fbexists = ["{}_{}".format(fluxbound.getReactionId(), fluxbound.getType()) for fluxbound in model.flux_bounds]
        ```

        """
        assert type(fluxbound) == FluxBound, '\nERROR: requires a FluxBound object, not something of type {}'.format(type(fluxbound))
        assert fluxbound.__objref__ is None, 'ERROR: object already bound to \"{}\", do you want to add a clone instead'.format(str(fluxbound.__objref__).split('to')[1][1:-1])
        if fluxbound.getId() in self.__global_id__:
            raise RuntimeError('Duplicate fluxbound ID detected: {}'.format(fluxbound.getId()))
        else:
            self.__pushGlobalId__(fluxbound.getId(), fluxbound)

        if fbexists is None:
            RL = self.getFluxBoundByReactionID(fluxbound.getReactionId(), 'lower')
            RU = self.getFluxBoundByReactionID(fluxbound.getReactionId(), 'upper')
            GO = True
            if RL != None and RL.is_bound == fluxbound.is_bound and RL.reaction == fluxbound.reaction:
                #print('Lower bound on reaction {} already exists'.format(fluxbound.getReactionId()))
                GO = False
            if RU != None and RU.is_bound == fluxbound.is_bound and RU.reaction == fluxbound.reaction:
                #print('Upper bound on reaction {} already exists'.format(fluxbound.getReactionId()))
                GO = False
            if __DEBUG__: print('Adding fluxbound: {}'.format(fluxbound.id))
            if GO:
                fluxbound.__objref__ = weakref.ref(self)
                self.flux_bounds.append(fluxbound)
            else:
                print('\"{}\" FluxBound for reaction {} exists, skipping'.format(fluxbound.is_bound, fluxbound.reaction))
        else:
            if "{}_{}".format(fluxbound.getReactionId(), fluxbound.getType()) not in fbexists:
                fluxbound.__objref__ = weakref.ref(self)
                self.flux_bounds.append(fluxbound)
            else:
                print('\"{}\" FluxBound for reaction {} exists, skipping'.format(fluxbound.is_bound, fluxbound.reaction))

    def addSpecies(self, species):
        """
        Add an instantiated Species object to the FBA model

        - *species* an instance of the Species class

        """
        assert isinstance(species, Species), '\nERROR: requires a Species object, not something of type {}'.format(type(species))
        assert species.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(species.__objref__).split('to')[1][1:-1])
        if __DEBUG__: print('Adding species: {}'.format(species.id))
        if species.getId() in self.__global_id__:
            raise RuntimeError('Duplicate species ID detected: {}'.format(species.getId()))
        else:
            self.__pushGlobalId__(species.getId(), species)
        species.__objref__ = weakref.ref(self)
        self.species.append(species)

    def addGene(self, gene):
        """
        Add an instantiated Gene object to the FBA model

        - *gene* an instance of the G class

        """
        # TODO: fix this whole gene thing, genes must use labels for gene names and id's for object search
        assert isinstance(gene, Gene), '\nERROR: requires a Gene object, not something of type {}'.format(type(gene))
        assert gene.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(gene.__objref__).split('to')[1][1:-1])
        if __DEBUG__: print('Adding Gene: {}'.format(gene.id))
        if gene.getId() in self.__global_id__:
            raise RuntimeError('Duplicate gene ID detected: {}'.format(gene.getId()))
        # Kill me
        elif gene.getId() in self.__genes_idx__:
            raise RuntimeError('Duplicate gene ID detected in index: {}'.format(gene.getId()))
        else:
            self.__pushGlobalId__(gene.getId(), gene)
        gene.__objref__ = weakref.ref(self)
        self.genes.append(gene)


    def addParameter(self, par):
        """
        Add an instantiated Parameter object to the model

        - *par* an instance of the Parameter class

        """
        assert isinstance(par, Parameter), '\nERROR: requires a Parameter object, not something of type {}'.format(type(par))
        assert par.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(par.__objref__).split('to')[1][1:-1])
        if __DEBUG__: print('Adding Parameter: {}'.format(par.id))
        if par.getId() in self.__global_id__:
            raise RuntimeError('Duplicate par ID detected: {}'.format(par.getId()))
        else:
            self.__pushGlobalId__(par.getId(), par)
        self.parameters.append(par)

    def addCompartment(self, comp):
        """
        Add an instantiated Compartment object to the CBM model

        - *comp* an instance of the Compartment class

        """
        assert isinstance(comp, Compartment), '\nERROR: requires a Compartment object, not something of type {}'.format(type(comp))
        assert comp.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(comp.__objref__).split('to')[1][1:-1])
        if __DEBUG__: print('Adding Compartment: {}'.format(comp.id))
        if comp.getId() in self.__global_id__:
            raise RuntimeError('Duplicate comp ID detected: {}'.format(comp.getId()))
        else:
            self.__pushGlobalId__(comp.getId(), comp)
        comp.__objref__ = weakref.ref(self)
        self.compartments.append(comp)

    def addReaction(self, reaction, create_default_bounds=False, silent=False):
        """
        Adds a reaction object to the model

        - *reaction* an instance of the Reaction class
        - *create_default_bounds* create default reaction bounds, irreversible 0 <= J <= INF, reversable -INF <= J <= INF
        - *silent* [default=False] if enabled this disables the printing of information messages

        """
        assert isinstance(reaction, Reaction), '\nERROR: requires a Reaction object, not something of type {}'.format(type(reaction))
        assert reaction.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(reaction.__objref__).split('to')[1][1:-1])

        if __DEBUG__: print('Adding reaction: {}'.format(reaction.id))
        if reaction.getId() in self.__global_id__:
            raise RuntimeError('Duplicate reaction ID detected: {}'.format(reaction.getId()))
        else:
            self.__pushGlobalId__(reaction.getId(), reaction)
        reaction.__objref__ = weakref.ref(self)
        for rr in reaction.reagents:
            rr.__objref__ = weakref.ref(self)
            if rr.getId() in self.__global_id__:
                raise RuntimeError('Duplicate reagent ID detected: {}'.format(reaction.getId()))
            else:
                self.__pushGlobalId__(rr.getId(), rr)
        self.reactions.append(reaction)
        if create_default_bounds:
            rid = reaction.getId()
            self.createReactionUpperBound(rid, numpy.inf)
            if reaction.reversible:
                self.createReactionLowerBound(rid, -numpy.inf)
                if not silent:
                    print('\nReaction \"{}\" bounds set to: -INF <= {} <= INF'.format(rid, rid))
            else:
                self.createReactionLowerBound(rid, 0)
                if not silent:
                    print('\nReaction \"{}\" bounds set to: 0 <= {} <= INF'.format(rid, rid))

    def addUserConstraint(self, pid, fluxes=None, operator='>=', rhs=0.0):
        """
        Add a user defined constraint to FBA model, this is additional to the automatically determined Stoichiometric constraints.

         - *pid* user constraint name/id, use `None` for auto-assign
         - *fluxes* a list of (coefficient, reaction id) pairs where coefficient is a float
         - *operator* is one of '=', '>=' or '<=' (< and > will be interpreted as >= or <=)
         - *rhs* a float

        """
        assert fluxes != None, '\nNo *fluxes* defined: a list of (coefficient, reaction id) pairs where coefficient is a float'
        if self.user_constraints == None:
            self.user_constraints = {}
        Rids = self.getReactionIds()
        fluxlist = []
        for J in fluxes:
            assert J[1] in Rids, '\n%s is not a valid reaction id' % J[1]
            assert J[1] not in fluxlist, '\nFluxes may only appear once per constraint:\n%s' % fluxes
            fluxlist.append(J[1])
        operator = operator.strip()
        assert operator in ['>', '<', '<=', '>=', '=', 'L', 'G', 'E'], '\n{} is not a valid operator'.format(operator)
        if operator == '=' or operator == 'E':
            operator = 'E'
        elif operator in ['>=', '>', 'G']:
            operator = 'G'
        elif operator in ['<=', '=<', '<', 'L']:
            operator = 'L'

        rhs = float(rhs)
        if pid is None:
            pid = 'uConstr%s' % (len(self.user_constraints)+1)
        self.user_constraints.update({pid : {'fluxes' : fluxes,
                                              'operator' : operator,
                                              'rhs' : rhs}
                                      })

    def deleteReactionAndBounds(self, rid):
        """
        Delete all reaction and bounds connected to reaction

         - *rid* a valid reaction id

        """

        Ridx = None
        Robj = None
        assert rid in self.getReactionIds(), '\nOh Hellooooooooooooo'
        for r in range(len(self.reactions)-1, -1, -1):
            if self.reactions[r].getId() == rid:
                Ridx = rid
                Robj = self.reactions.pop(r)
                self.__popGlobalId__(rid)
        Bounds = []
        for b in range(len(self.flux_bounds)-1, -1, -1):
            if self.flux_bounds[b].reaction == rid:
                Bounds.append(self.flux_bounds.pop(b))
                self.__popGlobalId__(Bounds[-1].getId())
        for re in Robj.reagents:
            S = self.getSpecies(re.getSpecies())
            if Ridx in S.reagent_of:
                S.reagent_of.remove(Ridx)
            self.__popGlobalId__(re.getId())
        Robj.reagents = []
        print('Deleting reaction {} and {} associated bounds'.format(Ridx, len(Bounds)))
        del Bounds, Robj
        #removed until I have a more secure way of doing this
        #self.__TRASH__.update({Ridx : {'react' : Robj.clone(), 'bnds' : [b.clone() for b in Bounds]}})


    def deleteObjective(self, objective_id):
        """
        Delete objective function:

         *objective_id* the id of the objective function. If objective_id is given  as 'active' then the active objective is deleted.

        """
        if objective_id == 'active':
            objective_id = self.getActiveObjective().getId()

        assert objective_id in self.getObjectiveIds(), '\nNo ....'
        for o in range(len(self.objectives)-1, -1, -1):
            if self.objectives[o].getId() == objective_id:
                Oobj = self.objectives.pop(o)
                for fo in Oobj.fluxObjectives:
                    self.__popGlobalId__(fo.getId())
                self.__popGlobalId__(objective_id)
        print('Deleting objective {}'.format(objective_id))
        del Oobj
        #removed until I have a more secure way of doing this
        #self.__TRASH__.update({objective_id : Oobj.clone()})

    #removed until I have a more secure way of doing this
    #def undeleteObjective(self, objective_id):
        #"""
        #Undeltes a deleted objective function:

         #- *objective_id* the id of an objeective function

        #"""

        #assert objective_id in self.__TRASH__, '\nNo deleted object of with this id'
        #self.addObjective(self.__TRASH__[objective_id])


    def deleteBoundsForReactionId(self, rid, lower=True, upper=True):
        """
        Delete bounds connected to reaction, rid

         - *rid* a valid reaction id
         - *upper* [default=True] delete the upper bound
         - *lower* [default=True] delete the lower bound

        """
        ##  Ridx = None
        ##  Robj = None
        Bounds = []
        for b in range(len(self.flux_bounds)-1, -1, -1):
            if self.flux_bounds[b].reaction == rid:
                delbound = False
                if upper and self.flux_bounds[b].getType() == 'upper':
                    delbound = True
                if lower and self.flux_bounds[b].getType() == 'lower':
                    delbound = True
                if lower and upper and self.flux_bounds[b].getType() == 'equality':
                    delbound = True
                if delbound:
                    fb = self.flux_bounds.pop(b)
                    self.__popGlobalId__(fb.getId())
                    Bounds.append(fb)
        print('Deleting {} bounds associated with reaction {}'.format(len(Bounds), rid))
        del Bounds

    #removed until I have a more secure way of doing this
    #def undeleteReactionAndBounds(self, rid):
        #"""
        #Undelete a reaction and bounds deleted with the **deleteReactionAndBounds** method

         #- *rid* a deleted reaction id

        #Please note this method is still experimental ;-)

        #"""
        #assert rid in self.__TRASH__, '\nOops I did it again ...'
        #Ridx = self.__TRASH__[rid]['react']
        #self.addReaction(Ridx, create_default_bounds=False)
        ## this is just while transitioning to new dynamic structures ... gone!
        ##for s_ in Ridx.getSpeciesIds():
            ##self.getSpecies(s_).setReagentOf(rid)
        #for b in self.__TRASH__[rid]['bnds']:
            #self.addFluxBound(b)
        #self.__TRASH__.pop(rid)

        #print('Undeleting reaction: {}'.format(rid))

    def createGeneAssociationsFromAnnotations(self, annotation_key='GENE ASSOCIATION', replace_existing=True):
        """
        Add genes to the model using the definitions stored in the annotation key. If this fails it tries some standard annotation
        keys: GENE ASSOCIATION, GENE_ASSOCIATION, gene_association, gene association.

         - *annotation_key* the annotation dictionary key that holds the gene association for the protein/enzyme
         - *replace_existing* [default=True] replace existing annotations, otherwise only new ones are added

        """
        # protein, assoc, id=None, name=None, gene_pattern='(\(\W*\w*\W*\))')
        if replace_existing:
            self.genes = []
            self.gpr = []
            self.__genes_idx__ = []
        gid = name = None
        g0 = len(self.genes)
        gpr0 = len(self.gpr)
        ga_keys = []
        for r_ in self.getReactionIds():
            GA = None
            R = self.getReaction(r_)
            ##  print r.annotation
            if annotation_key in R.annotation:
                GA = annotation_key
            elif 'GENE ASSOCIATION' in R.annotation:
                GA = 'GENE ASSOCIATION'
            elif 'GENE_ASSOCIATION' in R.annotation:
                GA = 'GENE_ASSOCIATION'
            elif 'gene_association' in R.annotation:
                GA = 'gene_association'
            elif 'gene association' in R.annotation:
                GA = 'gene association'
            if GA != None:
                self.createGeneProteinAssociation(r_, R.getAnnotation(GA), gid, name, update_idx=False)
                if GA not in ga_keys:
                    ga_keys.append(GA)

        print('INFO: used key(s) \'{}\''.format(ga_keys))
        self.__updateGeneIdx__()
        print('INFO: Added {} new genes and {} associations to model'.format(len(self.genes)- g0, len(self.gpr)- gpr0))


    def addGPRAssociation(self, gpr, update_idx=True):
        """
        Add a GeneProteinAssociation instance to the model

         - *gpr* an instantiated GeneProteinAssociation object

        """
        assert gpr.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(gpr.__objref__).split('to')[1][1:-1])

        if gpr.getId() in self.__global_id__:
            raise RuntimeError('Duplicate gpr ID detected: {}'.format(gpr.getId()))
        else:
            self.__pushGlobalId__(gpr.getId(), gpr)

        # will soon be terminated
        if update_idx:
            self.__updateGeneIdx__()

        gpr.__objref__ = weakref.ref(self)
        self.gpr.append(gpr)

    def createGeneProteinAssociation(self, protein, assoc, gid=None, name=None, gene_pattern=None, update_idx=True, altlabels=None):
        """
        Create and add a gene protein relationship to the model, note genes are mapped on protein objects which may or may not be reactions

         - *protein* in this case the reaction
         - *assoc* the COBRA style gene protein association
         - *gid* the unique id
         - *name* the optional name
         - *gene_pattern* deprecated, not needed anymore
         - *update_idx* update the model gene index, not used
         - *altlabels* [default=None] alternative labels for genes, default uses geneIds

        """
        if altlabels is None:
            altlabels = {}
        if assoc != '' and assoc is not None:
            if gid == None:
                gid = '{}_assoc'.format(protein)
            gpr = GeneProteinAssociation(gid, protein)
            self.addGPRAssociation(gpr)
            if name == None:
                name = gid
            gpr.setName(name)
            gpr.createAssociationAndGeneRefsFromString(assoc, altlabels)

    def createGeneProteinAssociationFromTree(self, protein, gprtree, gid=None, name=None, gene_pattern=None, update_idx=True, altlabels=None):
        """
        Create and add a gene protein relationship to the model, note genes are mapped on protein objects which may or may not be reactions

         - *protein* in this case the reaction
         - *gprtree* the CBMPy GPR dictionary tree
         - *gid* the unique id
         - *name* the optional name
         - *gene_pattern* deprecated, not needed anymore
         - *update_idx* update the model gene index, not used
         - *altlabels* [default=None] alternative labels for genes, default uses geneIds

        """
        if altlabels is None:
            altlabels = {}
        if gprtree != '' and gprtree is not None:
            if gid == None:
                gid = '%s_assoc' % protein
            gpr = GeneProteinAssociation(gid, protein)
            self.addGPRAssociation(gpr)
            if name == None:
                name = gid
            gpr.setName(name)
            gpr.createAssociationAndGeneRefsFromTree(gprtree, altlabels)

    def __updateGeneIdx__(self):
        """
        Update the gene name index

        """
        self.__genes_idx__ = [g.getId() for g in self.genes]

    def getAllProteinGeneAssociations(self, use_labels=False):
        """
        Returns a dictionary of the proteins associated with each gene

         - *use_labels* use V2 gene labels rather than ID's

        """
        prg = {}
        for gpr in self.gpr:
            for g in gpr.getGenes():
                if use_labels:
                    gid = g.getLabel()
                else:
                    gid = g.getId()
                if gid not in prg:
                    prg.update({gid : [gpr.protein]})
                else:
                    prg[gid].append(gpr.protein)
        return prg

    def getAllGeneProteinAssociations(self, use_labels=False):
        """
        Returns a dictionary of genes associated with each protein

         - *use_labels* use V2 gene labels rather than ID's

        """
        gprmap = {}
        for gpr in self.gpr:
            if gpr.protein not in gprmap:
                if use_labels:
                    gprmap.update({gpr.protein : gpr.getGeneLabels()})
                else:
                    gprmap.update({gpr.protein : gpr.getGeneIds()})
            else:
                if use_labels:
                    gprmap[gpr.protein].extend(gpr.getGeneLabels())
                else:
                    gprmap[gpr.protein].extend(gpr.getGeneIds())
        return gprmap

    def getGene(self, g_id):
        """
        Returns a gene object that has the identifier:

         - *gid* the gene identifier

        """
        out = None
        for g_ in self.genes:
            if g_.getId() == g_id:
                out = g_
                break
        return out

    def getGPRassociation(self, gpr_id):
        """
        Returns a gene protein association object that has the identifier:

         - *gpr_id* the gene protein identifier

        """
        out = None
        for g_ in self.gpr:
            if g_.getId() == gpr_id:
                out = g_
                break
        return out

    def getGPRforReaction(self, rid):
        """
        Return the GPR associated with the reaction id:

         - *rid* a reaction id

        """
        out = None
        assert self.getReaction(rid) is not None, '\nERROR: \"{}\" is not a valid reaction id'.format(rid)
        for gpr_ in self.gpr:
            if gpr_.getProtein() == rid:
                out = gpr_
                break
        return out

    def getGPRforReactionAsDict(self, rid, useweakref=True):
        """
        Return the GPR associated with the reaction id as a nested dictionary structure:

         - *rid* a reaction id

        """
        gprDict = {}
        try:
            createAssociationDictFromNode(ast.parse(self.getGPRforReaction(rid).getAssociationStr()).body[0], gprDict, self,\
                                      useweakref=useweakref, cntr=0)
        except SyntaxError:
            gprDict = {}
        return gprDict

    def getGPRIdAssociatedWithGeneId(self, gid):
        """
        Return the GPR(s) associated with the gene id:

         - *gid* a gene id

        """
        if self.getGene(gid) is None:
            print('INFO: invalid gene id: {}'.format(gid))
            return None
        out = []
        for gpr_ in self.gpr:
            if gid in gpr_.getGeneIds():
                out.append(gpr_.getId())
        if len(out) == 0:
            out = None
        return out

    def getGPRIdAssociatedWithGeneLabel(self, label):
        """
        Return the GPR Id's associated with the gene label:

         - *label* a gene label

        """
        gl = self.getGeneByLabel(label)
        if gl is None:
            print('INFO: invalid gene label: {}'.format(label))
            return None
        return self.getGPRIdAssociatedWithGeneId(gl.getId())

    def getReactionActivity(self, rid):
        """
        If there is a GPR and genes associated with the reaction ID then return either active=True or inactive=False
        Note if there is no gene associated information then this will return active.

         - *rid* a reaction id

        """
        out = True
        GPR = self.getGPRforReaction(rid)
        if GPR != None:
            out = GPR.isProteinActive()
        return out

    def getGeneIds(self, substring=None):
        """
        Returns a list of gene Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [g.getId() for g in self.genes]
        else:
            return [g.getId() for g in self.genes if substring in g.getId()]

    def getGeneObjects(self, substring=None):
        """
        Returns a list of gene objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [g for g in self.genes]
        else:
            return [g for g in self.genes if substring in g.getId()]

    def getGeneObjectsByLabel(self, substring=None):
        """
        Returns a list of gene objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the label

        """
        if substring == None:
            return [g for g in self.genes]
        else:
            return [g for g in self.genes if substring in g.getLabel()]

    def getGPRIds(self, substring=None):
        """
        Returns a list of GPR Id's, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [g.getId() for g in self.gpr]
        else:
            return [g.getId() for g in self.gpr if substring in g.getId()]

    def getGPRObjects(self, substring=None):
        """
        Returns a list of GPR objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [g for g in self.gpr]
        else:
            return [g for g in self.gpr if substring in g.getId()]

    def getGeneLabels(self, substring=None):
        """
        Returns a list of gene labels (locus tags), applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the label

        """
        if substring == None:
            return [g.getLabel() for g in self.genes]
        else:
            return [g.getLabel() for g in self.genes if substring in g.getLabel()]

    def getAllGeneActivities(self):
        """
        Returns a dictionary of genes (if defined) and whether they are active or not
        """
        out = {}
        for g_ in self.genes:
            out[g_.getId()] = g_.isActive()
        return out

    def getAllProteinActivities(self):
        """
        Returns a dictionary of reactions (if genes and GPR's are defined) and whether they are active or not

        """
        out = {}
        for g_ in self.gpr:
            out[g_.getProtein()] = g_.isProteinActive()
        return out

    def setAllProteinActivities(self, activites, lower=0.0, upper=0.0):
        """
        Given a dictionary of activities [rid : boolean] pairs set all the corresponding reactions:

         - *activities* a dictionary of [rid : boolean] pairs
         - *lower* [default=0.0] the lower bound of the deactivated flux
         - *upper* [default=0.0] the upper bound of the deactivated flux

        """
        for r_ in activites:
            R = self.getReaction(r_)
            if activites[r_] and not R.__is_active__:
                R.reactivateReaction()
            elif not activites[r_] and R.__is_active__:
                R.deactivateReaction(lower, upper)

    def setAllInactiveGPRBounds(self, lower=0.0, upper=0.0):
        """
        Set all reactions that are inactive (as determined by gene and gpr evaluation) to bounds:

         - *lower* [default=0.0] the new lower bound
         - *upper* [default=0.0] the new upper bound

        """
        gpract = self.getAllProteinActivities()
        for r_ in gpract:
            if not gpract[r_]:
                self.__gene_deactivated_reactions__[r_] = self.getReactionBounds(r_)
                self.setReactionBounds(r_, lower, upper)
                print('Setting reaction {} bounds to [{}, {}]'.format(r_, lower, upper))

    def resetAllInactiveGPRBounds(self):
        """
        Resets all reaction bounds modified by the ``cmod.setAllInactiveGeneReactionBounds()`` method to their previous values

        """
        for r_ in self.__gene_deactivated_reactions__:
            r = self.__gene_deactivated_reactions__.pop(r_)
            self.setReactionBounds(r_, r[1], r[2])
            print('Resetting reaction {} bounds to [{}, {}]'.format(r_, r[1], r[2]))

    def changeAllFluxBoundsWithValue(self, old, new):
        """
        Replaces all flux bounds with value "old" with a new value "new":

         - *old* value
         - *new* value

        """
        cntr = 0
        if new in ['inf', 'INF', 'INFINITY', '-inf', '-INF', '-INFINITY']:
            new = float(new)
        for b in range(len(self.flux_bounds)-1,-1,-1):
            if self.flux_bounds[b].value == old:
                ##  db = fba.flux_bounds.pop(b)
                self.flux_bounds[b].value = new
                cntr += 1
        print('\nReplaced {} \"{}\" bounds with \"{}\"'.format(cntr, old, new))

    def deleteAllFluxBoundsWithValue(self, value):
        """
        Delete all flux bounds which have a specified value:

         - *value* the value of the flux bound(s) to delete

        """
        cntr = 0
        for b in range(len(self.flux_bounds)-1,-1,-1):
            if self.flux_bounds[b].value == value:
                db = self.flux_bounds.pop(b)
                self.__popGlobalId__(db.getId())
                del db
                cntr += 1
        print('\nDeleted {} \"{}\" bounds'.format(cntr, value))

    def deleteCompartment(self, sid, check_components=True):
        """
        Deletes a compartment object with id. Returns True if the compartment is deleted, False if not. In addition if components were checked
        a list of id's that reference the compartment are also returned.

        - *sid* the compartment id
        - *check_components* [default=True] if  enabled check that no species or reactions makes
        use of the compartment, fail if it does.

        """

        if sid not in self.getCompartmentIds():
            print('\n\"{}\" is not a valid compartment id.\n'.format(sid))
            return False
        out = True
        if check_components:
            rc = [r.getId() for r in self.reactions if r.compartment == sid]
            sc = [s.getId() for s in self.species if s.compartment == sid]
            pc = [p.getId() for p in self.parameters if p.compartment == sid]
            mc = []

            if len(rc) > 0:
                err = ','.join(rc)
                print('Reactions: {} are located in compartment \"{}\"'.format(err, sid))
                out = False
            if len(sc) > 0:
                err = ','.join(sc)
                print('Species: {} are located in compartment \"{}\"'.format(err, sid))
                out = False
            if len(pc) > 0:
                err = ','.join(pc)
                print('Parameters: {} are located in compartment \"{}\"'.format(err, sid))
                out = False
            if sid == self.compartment:
                print('The model is located in compartment: {}'.format(err, sid))
                mc = [self.getId()]
                out = False

        if out:
            self.compartments.remove(self.getCompartment(sid))
            return True
        else:
            print('\nDeleteCompartment failed to delete compartment: {}\n'.format(sid))
            return out, sc+mc+pc+rc


    def deleteSpecies(self, sid, also_delete=None):
        """
        Deletes a species object with id

        - *sid* the species id
        - *also_delete* [default=None] only delete the species

        -- 'reagents' delete the species from the reactions it participates in as a **reagent**
        -- 'reactions' deletes the **reactions** that the species participates in

        """
        if also_delete is not None:
            rids = [r[0] for r in self.getFluxesAssociatedWithSpecies(sid)]
            if also_delete == 'reagents':
                for r in rids:
                    R = self.getReaction(r)
                    R.deleteReagentWithSpeciesRef(sid)
            elif also_delete == 'reactions':
                for r in rids:
                    self.deleteReactionAndBounds(r)
        SP = self.species.pop(self.getSpeciesIds().index(sid))
        self.__popGlobalId__(sid)
        del SP
        #removed until I have a more secure way of doing this
        #self.__TRASH__[sid] = SP.clone()

    #removed until I have a more secure way of doing this
    #def undeleteSpecies(self, sid):
        #"""
        #Undeltes a deleted species:

         #- *species* the species id

        #"""

        #assert sid in self.__TRASH__, '\nNo deleted object of with this id'
        #self.addSpecies(self.__TRASH__[sid])

    def deleteGene(self, gid, also_delete_gpr=True):
        """
        Deletes the gene object with gid. Note if you want to delete a gene by label (locus tag etc)
        use the deleteGeneByLabel() function.

        - *gid* the gene Id
        - *also_delete_gpr* [default=True] automatically delete GPR's that contain no gene references

        """
        print('DeleteGene is processing gene: {} ...'.format(gid))
        G = self.getGene(gid)
        if G is not None:
            assoc_gpr = self.getGPRIdAssociatedWithGeneId(gid)
            #print('DeleteGene associated GPRs: {}'.format(assoc_gpr))
            if assoc_gpr is not None:
                for gpr in assoc_gpr:
                    #print('DeleteGene is processing GPR: {}'.format(gpr))
                    GPR = self.getGPRassociation(gpr)
                    GPR.deleteGeneFromAssociation(gid)
                    if len(GPR.getGeneIds()) == 0 and also_delete_gpr:
                        self.deleteGPRAssociation(gpr)
            self.__popGlobalId__(gid)
            self.genes.pop(self.genes.index(G))
            #print('Gene: {}'.format(self.getGene(gid)))
        else:
            print('INFO: Gene Id \"{}\" does not exist'.format(gid))


    def deleteGeneByLabel(self, label, also_delete_gpr=True):
        """
        Deletes the gene object with label (b2003 etc).

        - *label* the gene with label to be deleted
        - *also_delete_gpr* [default=True] automatically delete GPR's that contain no gene references

        """

        gid = self.getGeneIdFromLabel(label)
        if gid is not None:
            self.deleteGene(gid, also_delete_gpr=also_delete_gpr)
        else:
            print('INFO: Gene label \"{}\" does not exist'.format(label))

    def deleteGPRAssociation(self, gprid):
        """
        Delete a GPR association with id

        - *gprid* the GPR association id

        """
        GPR = self.getGPRassociation(gprid)
        if GPR is not None:
            self.__popGlobalId__(gprid)
            self.gpr.pop(self.gpr.index(GPR))
        else:
            print('INFO: GPR Id \"{}\" does not exist'.format(gprid))

    def __pushGlobalId__(self, sid, obj):
        if not self.__ENABLE_GLOBAL_WEAKREF__:
            self.__global_id__[sid] = True
        else:
            self.__global_id__[sid] = obj

    def __popGlobalId__(self, sid):
        #if not self.__ENABLE_GLOBAL_WEAKREF__:
        x = self.__global_id__.pop(sid)
        del x

    def __changeGlobalId__(self, old, new, obj):
        if not self.__ENABLE_GLOBAL_WEAKREF__:
            self.__global_id__[new] = True
            self.__global_id__.pop(old)
        else:
            self.__global_id__[new] = obj
            self.__global_id__.pop(old)

    def deleteNonReactingSpecies(self, simulate=True):
        """
        Deletes all species that are not reagents (do not to take part in a reaction).
        *Warning* this deletion is permanent and greedy (not selective). Returns a list of (would be) deleted species

         - *simulate* [default=True] only return a list of the speciesId's that would have been deleted if False

        """
        active_reagents = []
        for R in self.reactions:
            for S in R.getSpeciesIds():
                if S not in active_reagents:
                    active_reagents.append(S)
            if __DEBUG__: print(active_reagents)
        if not simulate:
            print('Deleting non-reactive species', end=" ")
        deleted_species = []
        for S in range(len(self.species)-1,-1,-1):
            if self.species[S].getId() not in active_reagents:
                deleted_species.append(self.species[S].getId())
                if simulate:
                    pass
                else:
                    self.deleteSpecies(self.species[S].getId())
                if __DEBUG__: print(self.getSpeciesIds())
                print('.', end= " ")
        print(' ')
        return deleted_species

    def getCompartment(self, cid):
        """
        Returns a compartment object with *cid*

         - *cid* compartment ID

        """
        out = None
        for c in self.compartments:
            if c.getId() == cid:
                out = c
                break
        if self.compartments.count(cid) > 1:
            print('\nERROR: multiple compartments with id \"{}\" returning first'.format(cid))
        return out

    def getReaction(self, rid):
        """
        Returns a reaction object with *id*

         - *rid* reaction ID

        """
        out = None
        for r in self.reactions:
            if r.getId() == rid:
                out = r
                break
        return out

    def getSpecies(self, sid):
        """
        Returns a species object with *sid*

         - *sid* a specied ID

        """
        out = None
        for s in self.species:
            if s.getId() == sid:
                out = s
                break
        return out

    def getParameter(self, pid):
        """
        Returns a parameter object with pid

        """
        out = None
        for p in self.parameters:
            if p.getId() == pid:
                out = p
                break
        return out

    def getReactionBounds(self, rid):
        """
        Get the bounds of a reaction, returns a tuple of rid, lowerbound value, upperbound value and equality value (None means bound does not exist).

         - *rid* the reaction ID

        """
        lb = ub = eq = None
        lb = self.getFluxBoundByReactionID(rid, 'lower')
        ub = self.getFluxBoundByReactionID(rid, 'upper')
        eq = self.getFluxBoundByReactionID(rid, 'equality')
        if lb != None:
            if numpy.isinf(lb.value) or numpy.isreal(lb.value):
                lb = lb.value
            else:
                lb = float(lb)
        if ub != None:
            if numpy.isinf(ub.value) or numpy.isreal(ub.value):
                ub = ub.value
            else:
                ub = float(ub)
        if eq != None:
            if numpy.isinf(eq.value) or numpy.isreal(eq.value):
                eq = eq.value
            else:
                eq = float(eq)
        return rid, lb, ub, eq

    def getReactionLowerBound(self, rid):
        """
        Returns the lower bound of a reaction (it it exists) or None

         - *rid* the reaction ID

        """
        lb = eq = None
        lb = self.getFluxBoundByReactionID(rid, 'lower')
        if lb != None:
            if type(lb.value) != str and (numpy.isreal(lb.value) or numpy.isinf(lb.value)):
                lb = lb.value
            else:
                lb = float(lb.value)
        else:
            eq = self.getFluxBoundByReactionID(rid, 'equality')
            if eq != None:
                #print('\nINFO: Lower bound defined as an equality ({})'.format(rid))
                if numpy.isinf(eq.value) or numpy.isreal(eq.value):
                    lb = eq.value
                else:
                    lb = float(eq)
        return lb

    def getReactionUpperBound(self, rid):
        """
        Returns the upper bound of a reaction (it it exists) or None

         - *rid* the reaction ID

        """
        ub = eq = None
        ub = self.getFluxBoundByReactionID(rid, 'upper')
        if ub != None:
            if type(ub.value) != str and (numpy.isreal(ub.value) or numpy.isinf(ub.value)):
                ub = ub.value
            else:
                ub = float(ub.value)
        else:
            eq = self.getFluxBoundByReactionID(rid, 'equality')
            if eq != None:
                #print('\nINFO: Upper bound defined as an equality ({})'.format(rid))
                if numpy.isinf(eq.value) or numpy.isreal(eq.value):
                    ub = eq.value
                else:
                    ub = float(eq)
        return ub

    ##  def getBoundByName(self, rid, bound):
        ##  """
        ##  Return a FluxBound instance. Note this is an old name for the newer preferred method: `getFluxBoundByReactionID`

         ##  - *rid* the reaction ID
         ##  - *bound* the bound: 'upper', 'lower', 'equal'

        ##  """
        ##  print 'Deprecated: use *getFluxBoundByReactionID*'
        ##  return self.getFluxBoundByReactionID(rid, bound)

    def getFluxBoundByID(self, fid):
        """
        Returns a FluxBound with id

         - *fid* the fluxBound ID

        """
        c_ = None
        for c_ in self.flux_bounds:
            if c_.getId() == fid:
                return c_
        return c_

    def getFluxBoundByReactionID(self, rid, bound):
        """
        Returns a FluxBound instance

         - *rid* the reaction ID
         - *bound* the bound: 'upper', 'lower', 'equality'

        """
        c_ = None
        for c_ in self.flux_bounds:
            if c_.reaction == rid and c_.is_bound == bound:
                return c_

    def getFluxBoundsByReactionID(self, rid):
        """
        Returns all FluxBound instances connected to a reactionId as a tuple of valid
        (lower, upper, None) or (None, None, equality) or alternatively invalid (lower, upper, equality).

         - *rid* the reaction ID

         *under evaluation*

        """
        lower = None
        upper = None
        equality = None
        for c_ in self.flux_bounds:
            if c_.reaction == rid and c_.getType() == 'lower':
                lower = c_
            elif c_.reaction == rid and c_.getType() == 'upper':
                upper = c_
            elif c_.reaction == rid and c_.getType() == 'equality':
                equality = c_
            if lower != None and upper != None:
                return (lower, upper, None)
            elif equality != None:
                return (None, None, equality)
        print('\nSomething strange in function getFluxBoundsByReactionID, returns: {}{}{}'.format(lower, upper, equality))
        return (lower, upper, equality)

    def getCompartmentIds(self, substring=None):
        """
        Returns a list of compartment Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [c.getId() for c in self.compartments]
        else:
            return [c.getId() for c in self.compartments if substring in c.getId()]

    def getCompartmentObjects(self, substring=None):
        """
        Returns a list of compartment objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [c for c in self.compartments]
        else:
            return [c for c in self.compartments if substring in c.getId()]


    def getSpeciesIds(self, substring=None):
        """
        Returns a list of species Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [s.getId() for s in self.species]
        else:
            return [s.getId() for s in self.species if substring in s.getId()]

    def getSpeciesObjects(self, substring=None):
        """
        Returns a list of species objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring is None:
            return [s for s in self.species]
        else:
            return [s for s in self.species if substring in s.getId()]

    def getReactionIds(self, substring=None):
        """
        Returns a list of reaction Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [r.getId() for r in self.reactions]
        else:
            return [r.getId() for r in self.reactions if substring in r.getId()]

    def getReactionObjects(self, substring=None):
        """
        Returns a list of reaction objects, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [r for r in self.reactions]
        else:
            return [r for r in self.reactions if substring in r.getId()]

    def getReactionNames(self, substring=None):
        """
        Returns a list of reaction names, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the name

        """
        if substring == None:
            return [r.name for r in self.reactions]
        else:
            return [r.name for r in self.reactions if substring in r.name]

    def getFluxBoundIds(self, substring=None):
        """
        Returns a list of fluxbound Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [s.getId() for s in self.flux_bounds]
        else:
            return [s.getId() for s in self.flux_bounds if substring in s.getId()]

    def getObjectiveIds(self, substring=None):
        """
        Returns a list of objective function Ids, applies a substring search if substring is defined

         - *substring* search for this pattern anywhere in the id

        """
        if substring == None:
            return [s.getId() for s in self.objectives]
        else:
            return [s.getId() for s in self.objectives if substring in s.getId()]

    def getOptimalValue(self):
        """
        Returns the optimal value of the objective function

        """
        AO = self.getActiveObjective()
        print('Objective {}: \"{}\"'.format(AO.getId(), AO.operation))
        return AO.getValue()

    def getObjFuncValue(self):
        """
        Returns the objective function value

        """
        AO = self.getActiveObjective()
        print('Objective {}: \"{}\"'.format(AO.getId(), AO.operation))
        return AO.getValue()

    def getActiveObjective(self):
        """
        Returns the active objective object.

        """
        out = None
        try:
            out = self.objectives[self.activeObjIdx]
        except Exception:
            print('No active objective to get')
        return out

    def getActiveObjectiveStoichiometry(self):
        """
        Returns a list of (coefficient, flux_objective) tuples

        """
        out = None
        try:
            out = self.objectives[self.activeObjIdx].getFluxObjectiveData()
        except Exception:
            print('No active objective to get')
        return out


    def getActiveObjectiveReactionIds(self):
        """
        Returns the active objective flux objective reaction id's

        """
        out = None
        try:
            out = self.objectives[self.activeObjIdx].getFluxObjectiveReactions()
        except Exception:
            print('No active objective to get')
        return out

    def setActiveObjective(self, objId):
        idz = [o.id for o in self.objectives]
        if objId in idz:
            self.activeObjIdx = idz.index(objId)
            self.obj_func = self.objectives[self.activeObjIdx]
            print('Objective {} (idx={}) set active'.format(objId, self.activeObjIdx))
        else:
            print('Invalid ObjectiveId: {}'.format(objId))

    def setBoundValueByName(self, rid, value, bound):
        """
        Deprecated use setReactionBound

        Set a reaction bound
         - *rid* the reactions id
         - *value* the new value
         - *bound* this is either 'lower' or 'upper'

        """
        print("Deprecated: use *setReactionBound*")
        self.setReactionBound(rid, value, bound)

    def getGeneIdFromLabel(self, label):
        """
        Given a gene label it returns the corresponding Gene id or None

         - *label*

        """
        for g_ in self.genes:
            if g_.label == label:
                return g_.getId()
        return None

    def getGeneByLabel(self, label):
        """
        Given a gene label return the corresponding Gene object

         - *label*

        """
        for g_ in self.genes:
            if g_.label == label:
                return g_
        return None

    def setGeneInactive(self, g_id, update_reactions=False, lower=0.0, upper=0.0):
        """
        Effectively deletes a gene by setting it's inactive flag while optionally updating the GPR associated reactions

         - *g_id* a gene ID
         - *update_reactions* [default=False] update the associated reactions fluxbounds
         - *lower* [default=0.0] the deactivated reaction lower bound
         - *upper* [default=0.0] the deactivated reaction upper bound

        """
        g = self.getGene(g_id)
        if g != None:
            g.setInactive()
            if update_reactions:
                self.updateNetwork(lower, upper)
            #else:
                #self.__check_gene_activity__ = True
            return True
        else:
            return False

    def setGeneActive(self, g_id, update_reactions=False):
        """
        Effectively restores a gene by setting it's active flag

         - *g_id* a gene ID
         - *update_reactions* [default=False] update the associated reactions fluxbounds from the gene deletion bounds if they exist

        """
        g = self.getGene(g_id)
        if g != None:
            g.setActive()
            if update_reactions:
                self.updateNetwork()
            #else:
                #self.__check_gene_activity__ = True
            return True
        else:
            return False

    def updateNetwork(self, lower=0.0, upper=0.0, silent=False):
        """
        Update the reaction network based on gene activity. If reaction is deactivated then lower and upper bounds are used

         - *lower* [default=0.0] deactivated lower bound
         - *upper* [default=0.0] deactivated upper bound

        """
        for g_ in self.gpr:
            active = g_.isProteinActive()
            R = self.getReaction(g_.getProtein())
            if not active and R.__is_active__:
                R.deactivateReaction(lower, upper)
            elif active and not R.__is_active__:
                R.reactivateReaction()
        if not silent:
            print('Updating gene activity network ... done.')
        self.__check_gene_activity__ = False

    def resetAllGenes(self, update_reactions=False):
        """
        Resets all genes to their default activity state (normally on)

        - *update_reactions* [default=False] update the associated reactions fluxbounds from the gene deletion bounds if they exist

        """
        for g_ in self.genes:
            g_.resetActivity()
        self.updateNetwork()

    def setReactionBound(self, rid, value, bound):
        """
        Set a reaction bound

         - *rid* the reactions id
         - *value* the new value
         - *bound* this is either 'lower' or 'upper', or 'equal'

        """
        if rid not in self.getReactionIds():
            print('\nERROR setReactionBound: reaction id {} does not exist'.format(rid))
            return
        c2 = self.getFluxBoundByReactionID(rid, bound)
        # changed to no str() casting
        if c2 != None:
            c2.setValue(value)
            if __DEBUG__: print(c2.reaction, c2.operation, c2.value)
        else:
            # if the bound does not exist, create a new one
            if bound == 'upper':
                self.createReactionUpperBound(rid, value)
            elif bound == 'lower':
                self.createReactionLowerBound(rid, value)
            else:
                raise RuntimeError('\n%s is not a valid reaction name or \'%s\' bound does not exist ' % (rid,bound))

    def setReactionBounds(self, rid, lower, upper):
        """
        Set both the upper and lower bound of a reaction:

        - *rid* the good old reaction id
        - *lower* the lower bound value
        - *upper* the upper bound value

        """

        self.setReactionBound(rid, lower, 'lower')
        self.setReactionBound(rid, upper, 'upper')

    def setReactionLowerBound(self, rid, value):
        """
        Set a reactions lower bound (if it exists)

         - *rid* the reactions id
         - *value* the new value

        """
        self.setReactionBound(rid, value, 'lower')

    def setReactionUpperBound(self, rid, value):
        """
        Set a reactions upper bound (if it exists)

         - *rid* the reaction id
         - *value* the new value

        """
        self.setReactionBound(rid, value, 'upper')


    def getAllFluxBounds(self):
        """
        Returns a dictionary of all flux bounds [id:value]

        """
        out = {}
        for f_ in self.flux_bounds:
            fid = f_.getId()
            if fid in out:
                print('WARNING: duplicate fluxbound Id: {}'.format(fid))
            out[fid] = f_.getValue()
        return out

    def setAllFluxBounds(self, bounds):
        """
        DEPRECATED! use setFluxBoundsFromDict()

        Sets all the fluxbounds present in bounds

         - *bounds* a dictionary of [fluxbound_id : value] pairs (not per reaction!!!)

        """
        print('\nDEPRECATION WARNING: use setFluxBoundsFromDict instead of setAllFluxBounds\n')
        time.sleep(1)
        self.setFluxBoundsFromDict(bounds)

    def setFluxBoundsFromDict(self, bounds):
        """
        Sets all the fluxbounds present in bounds

         - *bounds* a dictionary of [fluxbound_id : value] pairs (not per reaction!!!)

        """

        fbids = self.getFluxBoundIds()
        for f_ in self.flux_bounds:
            if f_.getId() in fbids:
                f_.setValue(bounds[f_.getId()])


    def renameObjectIds(self, prefix=None, suffix=None, target='all', ignore=None):
        """
        This method is designed for target="all" other use may result in inconsistent models. Update: "species" and "reactions"
        should also work as advertised, please check results.

         - *prefix* [None] if supplied add as a prefix
         - *suffix* [None] if supplied add as a suffix
         - *target* ['all'] specify what class of objects to rename

          - 'species'
          - 'reactions'
          - 'bounds'
          - 'objectives'
          - 'all'

         - *ignore* [default=None] a list of id's to ignore

        """

        if ignore is None:
            ignore = []
        SUFFIX = PREFIX = False
        if suffix is None and prefix is None:
            raise RuntimeError("Either prefix or suffix must be defined")
        if prefix is not None:
            self.prefix = prefix
            PREFIX = True
        if suffix is not None:
            self.suffix = suffix
            SUFFIX = True

        ##  self.id = self.id+suffix
        if target == 'species' or target == 'all':
            for s in self.species:
                if s.getId() not in ignore:
                    if PREFIX:
                        #s.id = prefix+s.id
                        s.setId(prefix+s.getId())
                    if SUFFIX:
                        #s.id = s.id+suffix
                        s.setId(s.getId()+suffix)

        if target == 'reactions' or target == 'all':
            for s in self.reactions:
                if s.getId() not in ignore:
                    if PREFIX:
                        s.setId(prefix+s.getId())
                    if SUFFIX:
                        s.setId(s.getId()+suffix)
                #for r in s.reagents:
                    #if PREFIX:
                        #r.setId(prefix+r.getId())
                        #if r.species_ref not in ignore:
                            #r.species_ref = prefix+r.species_ref
                    #if SUFFIX:
                        #r.setId(r.getId()+suffix)
                        #if r.species_ref not in ignore:
                            #r.species_ref = r.species_ref+suffix

        if target == 'bounds' or target == 'all':
            for s in self.flux_bounds:
                if PREFIX:
                    s.setId(prefix+s.getId())
                if SUFFIX:
                    s.setId(s.getId()+suffix)

        if target == 'objectives' or target == 'all':
            for s in self.objectives:
                if PREFIX:
                    s.setId(prefix+s.getId())
                if SUFFIX:
                    s.setId(s.getId()+suffix)
                for f in s.fluxObjectives:
                    if PREFIX:
                        f.setId(prefix+f.getId())
                        #if f.reaction not in ignore:
                            #f.reaction = prefix+f.reaction
                    if SUFFIX:
                        f.setId(f.getId()+suffix)
                        #if f.reaction not in ignore:
                            #f.reaction = f.reaction+suffix
        self.buildStoichMatrix()

    def setPrefix(self, prefix, target):
        """
        This is alpha stuff, target can be:

         - 'species'
         - 'reactions'
         - 'constraints'
         - 'objectives'
         - 'all'

        """

        print('\nDEPRECATION WARNING: setPrefix will be removed, please use the new cmod.renameObjectIds() method instead\n')
        time.sleep(1.0)

        if target == 'constraints':
            target = 'bounds'
        self.renameObjectIds(prefix=prefix, suffix=None, target=target, ignore=[])


    def setSuffix(self, suffix, target):
        """
        This is alpha stuff, target can be:

         - 'species'
         - 'reactions'
         - 'constraints'
         - 'objectives'
         - 'all'

        """

        print('\nDEPRECATION WARNING: setSuffix will be removed, please use the new cmod.renameObjectIds() method instead\n')
        time.sleep(1.0)

        if target == 'constraints':
            target = 'bounds'
        self.renameObjectIds(prefix=None, suffix=suffix, target=target, ignore=[])


    def setObjectiveFlux(self, rid, coefficient=1, osense='maximize', delete_objflx=True):
        """
        Set single target reaction flux for the current active objective function.

         - *rid* a string containing a reaction id
         - *coefficient* [default=1] an objective flux coefficient
         - *osense* the optimization sense must be **maximize** or **minimize**
         - *delete_objflx* [default=True] delete all existing fluxObjectives in the active objective function
        """
        assert rid in self.getReactionIds(), '\n%s is not a reaction\n%s' % (rid, self.getReactionIds())
        osense = osense.lower()
        if osense == 'max': osense = 'maximize'
        if osense == 'min': osense = 'minimize'
        if osense in ['maximise', 'minimise']:
            osense = osense.replace('se', 'ze')
        assert osense in ['maximize', 'minimize'], "\nosense must be ['maximize', 'minimize'] not %s" % osense
        if delete_objflx:
            self.objectives[self.activeObjIdx].deleteAllFluxObjectives()
        FO = FluxObjective('{}_{}_fluxobj'.format(self.objectives[self.activeObjIdx].getId(), rid), rid, coefficient)
        self.objectives[self.activeObjIdx].addFluxObjective(FO)
        self.objectives[self.activeObjIdx].operation = osense

    def sortReactionsById(self):
        """
        Sorts the reactions by Reaction.id uses the python string sort

        """

        reactions2 = []
        reactId = self.getReactionIds()
        reactId.sort()
        for r in reactId:
            reactions2.append(self.getReaction(r))
        self.reactions = reactions2
        self.__updateGeneIdx__()
        for gpr_ in self.gpr:
            gpr_.buildEvalFunc()
        del reactions2

    def sortSpeciesById(self):
        """
        Sorts the reaction list by Reaction.id uses the python string sort

        """
        raise NotImplementedError("\nCall developer about this")

    def getExchangeReactions(self):
        """
        Returns reaction instances where the 'is_exchange' attribute set to True. This is by default
        reactions that contain a boundary species.

        """
        return [r for r in self.reactions if r.is_exchange]

    def getExchangeReactionIds(self):
        """
        Returns id's of reactions where the 'is_exchange' attribute set to True. This is by default
        reactions that contain a boundary species.

        """
        return [r.getId() for r in self.reactions if r.is_exchange]

    def getReactionValues(self, only_exchange=False):
        """
        Returns a dictionary of ReactionID : ReactionValue pairs:

         - *only_exchange* [default=False] only return the reactions labelled as exchange

        """
        out = {}
        if only_exchange:
            exrids = set(self.getExchangeReactionIds())
        for r in self.reactions:
            if only_exchange and not r.getId() in exrids:
                continue
            out[r.getId()] = r.getValue()
        return out

    def getSolutionVector(self, names=False):
        """
        Return a vector of solution values

         - *names* [default=False] if True return a solution vector and list of names

        """
        self.buildStoichMatrix()
        J = numpy.array([self.getReaction(r).getValue() for r in self.N.col])
        if not names:
            return J
        else:
            return J, tuple(self.N.col)

    def getReversibleReactionIds(self):
        """
        Return a list of reversible reaction Id's

        """
        output = []
        for r_ in self.reactions:
            if r_.reversible:
                output.append(r_.getId())
        return output

    def getIrreversibleReactionIds(self):
        """
        Return a list of irreversible reaction Id's

        """
        output = []
        for r_ in self.reactions:
            if not r_.reversible:
                output.append(r_.getId())
        return output

    def getBoundarySpeciesIds(self, rid=None):
        """
        Return all boundary species associated with reaction

         - rid [default=None] by default return all boundary species in a model, alternatively a string containing a reaction id or list of reaction id's

        """
        if rid is None:
            rid = self.getReactionIds()
        elif type(rid) == list:
            pass
        else:
            rid = [rid]
        bspec = []
        for r_ in rid:
            R = self.getReaction(r_)
            for s_ in R.getSpeciesIds():
                if s_ not in bspec and self.getSpecies(s_).is_boundary:
                    bspec.append(s_)
        return bspec

    def findFluxesForConnectedSpecies(self, metab):
        """
        Returns a list of (reaction, flux value) pairs that this metabolite appears as a reagent of

         - *metab* the metabolite name

        """
        print('DEPRECATION WARNING: please use cmod.getFluxesAssociatedWithSpecies()')
        return self.getFluxesAssociatedWithSpecies(metab)

    def getFluxesAssociatedWithSpecies(self, metab):
        """
        Returns a list of (reaction, flux value) pairs that this metabolite appears as a reagent in

         - *metab* the metabolite name

        """
        output = []
        assert self.getSpecies(metab) is not None, '\n%s is not a valid metabolite name!' % metab
        for r_ in self.getSpecies(metab).isReagentOf():
            # if cmod.getReaction(r_).getValue() != 0.0:
            output.append([r_, self.getReaction(r_).getValue()])
            #print r_, cmod.getReaction(r_).getValue()
        return output

    def getReactionIdsAssociatedWithSpecies(self, metab):
        """
        Returns a list of (reaction, flux value) pairs that this metabolite appears as a reagent in

         - *metab* the metabolite name

        """
        output = []
        assert self.getSpecies(metab) != None, '\n%s is not a valid metabolite name!' % metab
        for r_ in self.getSpecies(metab).isReagentOf():
            # if cmod.getReaction(r_).getValue() != 0.0:
            output.append([r_, self.getReaction(r_).getValue()])
            #print r_, cmod.getReaction(r_).getValue()
        return output


    def getFluxesAssociatedWithCompartments(self, compartments):

        """
        Determines all reactions and flux values associated with a list of
        compartments. This function can be used to find all transport reactions
        between compartments, e.g. the cytosol and mitochondria. If the
        compartment IDs are 'cyt' and 'mit', respectively, you can call
        "your_model.getFluxesAssociatedWithCompartments(['cyt', 'mit'])"
        to get all fluxes between these compartments.

        *compartments*: a list or set of compartment IDs.
                        To check the existing compartment IDs in your model
                        call "your_model.getCompartmentIds()"

        :returns a dictionary with reaction IDs as keys and corresponding
        flux values as values
        """

        # check whether provided compartment ID's are valid
        if (not set(compartments).issubset(self.getCompartmentIds()) or
            not isinstance(compartments, (list, set))):
            raise ValueError("Please provide valid compartment IDs as a list or"
                             " set!")

        compartments = set(compartments)

        # check whether provided compartments are identical with the ones of
        # the reagents of a reaction; if so, add it to dictionary along with
        # the flux value
        return {ri: self.getReaction(ri).getValue() for ri in
                self.getReactionIds() if compartments ==
                set(si.getCompartmentId() for si in
                self.getReaction(ri).getSpeciesObj())}


    def splitEqualityFluxBounds(self):
        """
        Splits any equalit flux bounds into lower and upper bounds.

        """
        ebs = []
        for b_ in self.flux_bounds:
            if b_.getType() == 'equality':
                ebs.append(b_.getId())
        for e_ in ebs:
            old = self.getFluxBoundByID(e_)
            R = old.getReactionId()
            V = old.getValue()
            self.deleteBoundsForReactionId(R)
            self.createReactionLowerBound(R, V)
            self.createReactionUpperBound(R, V)
        del ebs

    def buildStoichMatrix(self, matrix_type='numpy', only_return=False):
        """
        Build the stoichiometric matrix N and additional constraint matrix CN (if required)

         - *matrix_type* [default='numpy'] the type of matrix to use to generate constraints

           - *numpy* a NumPy matrix default
           - *sympy* a SymPy symbolic matrix, if available note the denominator limit can be set in ``CBModel.__CBCONFIG__['SYMPY_DENOM_LIMIT'] = 10**12``
           - *scipy_csr* create using NumPy but store as SciPy csr_sparse

        - *only_return* [default=False] **IMPORTANT** only returns the stoichiometric matrix and constraint matrix (if required),
          does not update the model

        """
        if __DEBUG__: print('Species:', self.getSpeciesIds())
        var_spec = [s for s in self.species if not s.is_boundary]
        var_spec_id = [s.getId() for s in self.species if not s.is_boundary]
        reac_id = self.getReactionIds()
        if len(var_spec_id) != len(set(var_spec_id)):
            print('\nBUILD STOICHIOMETRY WARNING: duplicate species IDs detected!\n')
        if len(reac_id) != len(set(reac_id)):
            print('\nBUILD STOICHIOMETRY WARNING: duplicate reaction IDs detected!\n')

        num_col = len(self.reactions)
        num_row = len(var_spec)

        sym_dlim = __CBCONFIG__['SYMPY_DENOM_LIMIT']
        if __DEBUG__: print('N-dimension = (%s, %s)' % (num_row, num_col))
        SYMGO = False
        SCIGO = False
        if matrix_type == 'scipy_csr':
            if HAVE_SCIPY:
                spsrc = {}
                spcol = []
                sprow = []
                spdata = []
                RHS = numpy.zeros(num_row)
                SCIGO = True
                print('\nINFO: using SciPy sparse for N')
            else:
                raise RuntimeError('\nSciPy required for scipy_csr matrices')
        elif matrix_type == 'sympy':
            if HAVE_SYMPY:
                N = sympy.zeros(num_row, num_col)
                RHS = [sympy.numbers.Zero() for i in range(num_row)]
                SYMGO = True
                print('\nINFO: using SymPy symbolic for N')
            else:
                raise RuntimeError('\nSymPy required for SymPy matrices')
        else:
            N = numpy.zeros((num_row, num_col))
            RHS = numpy.zeros(num_row)

        for c in range(num_col):
            if __DEBUG__: print(self.reactions[c].getId())
            if __DEBUG__: print(self.reactions[c].getStoichiometry())
            for reag in self.reactions[c].getStoichiometry():
                if reag[1] in var_spec_id:
                    if __DEBUG__: print('{}: setting reagent {} to {} (idx={},{})'.format(self.reactions[c].getId(), reag[1], reag[0], var_spec_id.index(reag[1]), c))
                    r = var_spec_id.index(reag[1])
                    if SYMGO:
                        if N[r, c] == 0.0:
                            N[r, c] = sympy.Rational(reag[0]).limit_denominator(sym_dlim)
                        else:
                            N[r, c] = N[r, c] + sympy.Rational(reag[0]).limit_denominator(sym_dlim)
                    elif SCIGO:
                        # we can consider using a threshold here, but let's be strict for now
                        if (r, c) in spsrc:
                            if spsrc[(r, c)] == 0.0 and reag[0] == 0.0:
                                pass
                            else:
                                spsrc[(r, c)] = spsrc[(r, c)] + reag[0]
                        else:
                            spsrc[(r, c)] = reag[0]
                    elif N[r, c] == 0.0:
                        N[r, c] = reag[0]
                    else:
                        N[r, c] = N[r, c] + reag[0]
        if SCIGO:
            rck = list(spsrc.keys())
            rck.sort()
            for r, c in rck:
                sprow.append(r)
                spcol.append(c)
                spdata.append(spsrc[(r,c)])
            N = csr_matrix((spdata, (sprow, spcol)), shape=(num_row, num_col), dtype='d')
            N = StructMatrixLP(N, list(range(num_row)), list(range(num_col)), row=var_spec_id, col=reac_id, rhs=RHS)
            N.array.eliminate_zeros()
        else:
            N = StructMatrixLP(N, list(range(num_row)), list(range(num_col)), row=var_spec_id, col=reac_id, rhs=RHS)

        # build and append additional constraint matric
        CM = None
        if self.user_constraints != None:
            crows = list(self.user_constraints)
            crows.sort()
            ccols = reac_id
            cnum_col = len(ccols)
            cnum_row = len(crows)
            Coperators = ['E']*cnum_row

            if SCIGO:
                spcol = []
                sprow = []
                spdata = []
                CRHS = numpy.zeros(cnum_row)
                print('\nINFO: using SciPy sparse for CM')
            elif SYMGO:
                CM = sympy.zeros(cnum_row, cnum_col)
                CRHS = [sympy.numbers.Zero() for i in range(cnum_row)]
                print('INFO: using SymPy for CM')
            else:
                CM = numpy.zeros((cnum_row, cnum_col))
                CRHS = numpy.zeros(cnum_row)
            for cs in range(cnum_row):
                for flx in self.user_constraints[crows[cs]]['fluxes']:
                    tcol = ccols.index(flx[1])
                    if SYMGO:
                        CM[cs, tcol] = sympy.Rational(flx[0]).limit_denominator(sym_dlim)
                    elif SCIGO:
                        # we can consider using a threshold here, but let's be strict for now
                        if float(flx[0]) == 0.0:
                            pass
                        else:
                            sprow.append(cs)
                            spcol.append(tcol)
                            spdata.append(float(flx[0]))
                    else:
                        CM[cs, tcol] = float(flx[0])

                    Coperators[cs] = self.user_constraints[crows[cs]]['operator']
                    if SYMGO:
                        CRHS[cs] = sympy.Rational(self.user_constraints[crows[cs]]['rhs']).limit_denominator(sym_dlim)
                    else:
                        CRHS[cs] = float(self.user_constraints[crows[cs]]['rhs'])

            if matrix_type == 'scipy_csr':
                CM = csr_matrix((spdata, (sprow, spcol)), shape=(cnum_row, cnum_col), dtype='d')
                CM = StructMatrixLP(CM, list(range(cnum_row)), list(range(cnum_col)), row=crows, col=ccols,\
                                    rhs=CRHS, operators=Coperators)
                CM.array.eliminate_zeros()
                spcol = []
                sprow = []
                spdata = []
            else:
                CM = StructMatrixLP(CM, list(range(cnum_row)), list(range(cnum_col)), row=crows, col=ccols,\
                            rhs=CRHS, operators=Coperators)
        if not only_return:
            self.N = N
            if CM != None:
                self.CM = CM
        else:
            print('INFO: Model stoichiometry not updated.')
            if CM != None:
                return (N, CM)
            else:
                return N


    def createSingleGeneEffectMap(self):
        """
        This takes a model and analyses the logical gene expression patterns. This only needs to be done once,
        the result is a dictionary that has boolean effect patterns as keys and the (list of) genes that give rise to
        those patterns as values. This map is used by the single gene deletion method for further analysis.

        Note this dictionary can also be stored and retrieved separately as long as the model structure is not changed i.e.
        the gene associations themselves or order of reactions (stored as the special entry 'keyJ').

        Stored as self.__single_gene_effect_map__

        """
        fba2 = self.clone()
        #cbm.analyzeModel(fba2)

        geneIds = fba2.getGeneIds()
        wtpr = fba2.getAllProteinActivities()
        wtrn = fba2.getReactionIds()
        wtrn.sort()

        cpress = {binHash(wtrn, wtpr) : ['wt']}

        #print(cpress)

        geneidcache = {}

        t0a = time.time()
        DEBUG = False

        print(time.strftime('%H:%M:%S'))
        cntr = 1
        for g_ in geneIds:
            if g_ not in geneidcache:
                geneidcache[g_] = fba2.getGene(g_)
            geneidcache[g_].setInactive()
            new = fba2.getAllProteinActivities()
            bhash = binHash(wtrn, new)
            if bhash in cpress:
                cpress[bhash].append(g_)
            else:
                cpress[bhash] = [g_]
            geneidcache[g_].setActive()
            if time.time() - t0a > 20.0:
                print(cntr, time.strftime('%H:%M:%S'))
                t0a = time.time()
            cntr += 1
        print(cntr, time.strftime('%H:%M:%S'))
        cpress['keyJ'] = tuple(wtrn)
        self.__single_gene_effect_map__ = cpress
        del geneidcache, fba2, cpress

    def testGeneProteinAssociations(self):
        """
        This method will test the GeneProtein associations and return a list of protein, association pairs

        """
        out = []
        for g_ in self.gpr:
            try:
                g_.evalAssociation()
            except Exception:
                out.append((g_.protein, g_.assoc))
                print('Error in GPR for reaction: {}'.format(g_.protein))
        return out

    def exportFVAdata(self):
        """
        Export the fva data as an array and list of reaction id's

        """

        # investigate the use of numpy index arrays here
        raise NotImplementedError

    def importUserConstraints(self, filename):
        """
        Exports user constraints in json

        """
        if filename is None:
            return False
        F = open(filename, 'r')
        #din = json.load(F)
        self.user_constraints = json.load(F)
        #key = os.path.split(self.sourcefile)[-1]
        #F.close()
        #if key in din:
            #self.user_constraints = din[key]
        #else:
            #print('ERROR: constraints refers to file \"{}\" whereas this is file \"{}\"'.format(list(din.keys())[0], key))
            #return False
        return True

    def exportUserConstraints(self, filename):
        """
        Exports user constraints in json

        """
        F = open(filename, 'w')
        #json.dump({os.path.split(self.sourcefile)[-1] : self.user_constraints}, F, indent=2)
        json.dump(self.user_constraints, F, indent=2)
        F.close()

    def createGroup(self, gid):
        """
        Create an empty group with

         - *gid* the unique group id

        """
        self.addGroup(Group(gid))

    def addGroup(self, obj):
        """
        Add an instantiated group object to the model

         - *obj* the Group instance

        """
        assert obj.__objref__ is None, 'ERROR: object already bound to \"{}\", add a clone instead'.format(str(obj.__objref__).split('to')[1][1:-1])
        if obj.getId() not in self.getGroupIds():
            self.groups.append(obj)
        else:
            print('ERROR: Group with id \"{}\" already exists.'.format(obj.getId()))
            del obj

    def deleteGroup(self, gid):
        """
        Delete a group with

         - *gid* the unique group id

        """
        gids = self.getGroupIds()
        if gid in gids:
            self.groups.pop(gids.index(gid))
        else:
            print('ERROR: Group with id \"{}\" does not exist.'.format(gid))

    def getGroup(self, gid):
        """
        Return a group with

         - *gid* the unique group id

        """
        for g_ in self.groups:
            if g_.getId() == gid:
                return g_
        return None

    def getGroupIds(self):
        """
        Get all group ids

        """
        return [g.getId() for g in self.groups]

    def getGroupNames(self):
        """
        Get all group names

        """
        return [g.getName() for g in self.groups]

    def getGroupMembership(self):
        """
        Returns group membership of items in groups. Returns {object_id: ['group_id1', 'group_id2']}

        """
        grps = {}
        for g in self.groups:
            gid = g.getId()
            for mid in g.getMemberIDs():
                if mid in grps:
                    grps[mid].append(gid)
                else:
                    grps[mid] = [gid]
        return grps

    def emptyUndelete(self):
        """
        Empties the undelete cache

        """
        self.__TRASH__.clear()


class Objective(Fbase):
    """
    An objective function

    """
    fluxObjectives = None
    ##  fluxObjectiveNames = None
    operation = None
    value = None
    solution = None

    def __init__(self, pid, operation):
        pid = str(pid)
        self.setId(pid)

        if operation.lower() in ['maximize', 'maximise', 'max']:
            self.operation = 'maximize'
        elif operation.lower() in ['minimize', 'minimise', 'min']:
            self.operation = 'minimize'
        else:
            print('WARNING: Invalid operation: {}'.format(operation))
        self.fluxObjectives = []
        self.compartment = None
        self.__delattr__('compartment')

    def setOperation(self, operation):
        """
        Sets the objective operation (sense)

         - *operation* [default='maximize'] one of 'maximize', 'maximise', 'max', 'minimize', 'minimise', 'min'

        """
        if operation.lower() in ['maximize', 'maximise', 'max']:
            self.operation = 'maximize'
        elif operation.lower() in ['minimize', 'minimise', 'min']:
            self.operation = 'minimize'
        else:
            print('WARNING: Invalid operation: {}'.format(operation))

    def getOperation(self):
        """
        Returns the operation or sense of the objective
        """
        return self.operation

    def addFluxObjective(self, fobj, override=False):
        """
        Adds a FluxObjective instance to the Objective

         - *fobj* the FluxObjective object
         - *override* [default=False] override pushing the global id map, this should never be used

        """
        if fobj.getId() in self.getFluxObjectiveIDs():
            print('\nWARNING: a flux objective with id \"{}\" already exists ... not adding!\n'.format(fobj.getId()))
            return
        if not override:
            self.__objref__().__pushGlobalId__(fobj.getId(), fobj)
        self.fluxObjectives.append(fobj)

    def createFluxObjectives(self, fluxlist):
        """
        Create and add flux objective objects to this objective function.

         - *fluxlist* a list of one or more ('coefficient', 'rid') pairs

        """
        FOreact = self.getFluxObjectiveReactions()
        for J in fluxlist:
            if J[1] not in FOreact:
                fid = '{}_{}_fobj'.format(self.getId(), J[1])
                self.addFluxObjective(FluxObjective(fid, J[1], J[0]))
            else:
                print('\nObjective {} already contains flux {} ... skipping!\n'.format(self.getId(), J[1]))

    def deleteAllFluxObjectives(self):
        """
        Delete all flux objectives

        """
        for fo in self.fluxObjectives:
            self.__objref__().__popGlobalId__(fo.getId())
        self.fluxObjectives = []

    def getFluxObjectiveIDs(self):
        """
        Returns a list of ObjectiveFlux ids, for the reaction id's use *getFluxObjectiveReactions()*
        or for coefficient, fluxobjective pairs use *getFluxObjectiveData()*

        """
        return [f.getId() for f in self.fluxObjectives]

    def getFluxObjectiveForReaction(self, rid):
        """
        Returns the FluxObjective associated with the suplied rid. If there is more than fluxObjective associated with a reaction (illegal)
        then a list of fluxObjectives is returned.

         *rid* a reaction id

        """
        fo = None
        for fo_ in self.fluxObjectives:
            if fo_.reaction == rid:
                if fo == None:
                    fo = fo_
                elif type(fo) == list:
                    fo.append(fo_)
                    print('\nWARNING: multiple fluxObjectives match rid: {}\n'.format(rid))
                else:
                    fo = [fo]
                    fo.append(fo_)
                    print('\nWARNING: multiple fluxObjectives match rid: {}\n'.format(rid))
        return fo

    def getFluxObjectiveReactions(self):
        """
        Returns a list of reactions that are used as FluxObjectives

        """
        return [f.reaction for f in self.fluxObjectives]

    def getFluxObjectiveData(self):
        """
        Returns a list of ObjectiveFunction components as (coefficient, flux) pairs

        """
        return [(f.coefficient, f.reaction) for f in self.fluxObjectives]

    def getFluxObjective(self, foid):
        """
        Return the flux objective with id.

         - *foid* the flux objective id returns either an object or a list if there are multiply defined flux objectives

        """
        fo = None
        for fo_ in self.fluxObjectives:
            if fo_.getId() == foid:
                if fo == None:
                    fo = fo_
                elif type(fo) == list:
                    fo.append(fo_)
                    print('ERROR: multiple FluxBounds have id: {}'.format(foid))
                else:
                    fo = [fo]
                    fo.append(fo_)
                    print('ERROR: multiple FluxBounds have id: {}'.format(foid))
        return fo


    def getFluxObjectives(self):
        """
        Returns the list of FluxObjective objects.

        """
        return self.fluxObjectives

    def getValue(self):
        """
        Returns the current value of the attribute (input/solution)
        """
        return self.value

    def setValue(self, value):
        """
        Sets the attribute ''value''
        """
        self.value = value

class FluxObjective(Fbase):
    """
    A weighted flux that appears in an objective function

    NOTE: reaction is a string containing a reaction id
    """
    reaction = None
    coefficient = None

    def __init__(self, pid, reaction, coefficient=1):
        pid = str(pid)
        self.setId(pid)

        self.reaction = reaction
        self.coefficient = coefficient
        self.annotation = {}
        self.compartment = None
        self.__delattr__('compartment')

    def getReactionId(self):
        return self.reaction

    def getCoefficient(self):
        return self.coefficient

    def setReactionId(self, reaction):
        self.reaction = reaction

    def setCoefficient(self, coefficient):
        self.coefficient = coefficient

class Compartment(Fbase):
    """A compartment"""
    size = None
    dimensions = None
    volume = None

    def __init__(self, pid, name=None, size=1, dimensions=3, volume=None):
        """
         - *cid* comapartment id
         - *name* optional compartment name
         - *size* compartment size, [default=1]
         - *dimensions* compartment spatial dimensions [default=3]
        """
        pid = str(pid)
        if pid == self.id:
            return

        if not self.__checkId__(pid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(pid))

        self.id = pid

        if name != None:
            self.setName(name)
        else:
            self.setName(pid)
        self.size = size
        if volume == None:
            self.volume = size
        self.dimensions = dimensions
        self.annotation = {}


    def setId(self, fid):
        """
        Sets the object Id

         - *fid* a valid c variable style id string

         Reimplements @FBase.setId()

        """

        fid = str(fid)
        if fid == self.id:
            return

        if not self.__checkId__(fid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(fid))

        if self.__objref__ is not None:
            if fid not in self.__objref__().__global_id__:
                self.__objref__().__changeGlobalId__(self.id, fid, self)
                oldid = self.id
                self.id = fid
                for s in self.__objref__().species:
                    if s.compartment == oldid:
                        s.setCompartmentId(fid)
                for r in self.__objref__().reactions:
                    if r.compartment == oldid:
                        r.setCompartmentId(fid)
            else:
                print('ERROR: setId() - object with id \"{}\" already exists ... ID *not* set.'.format(fid))
        else:
            self.id = fid

    def containsSpecies(self):
        """
        Lists the species contained in this compartment

        """
        out = []
        if self.__objref__ != None:
            out = [s.getId() for s in self.__objref__().species if s.compartment == self.getId()]
        return out

    def containsReactions(self):
        """
        Lists the species contained in this compartment

        """
        out = []
        if self.__objref__ != None:
            out = [r.getId() for r in self.__objref__().reactions if r.compartment == self.getId()]
        return out

    def getSize(self):
        """
        Get the compartment size

        """
        return self.size

    def getDimensions(self):
        """
        Get the compartment dimensions

        """
        return self.dimensions

    def setSize(self, size):
        """
        Set the compartment size

         - *size* the new compartment size

        """
        self.size = size

    def setDimensions(self, dimensions):
        """
        Get the compartment dimensions

         - *dimensions* set the new compartment dimensions

        """
        self.dimensions = dimensions

class GroupMemberAttributes(Fbase):
    """
    Contains the shared attributes of the group members (equivalent to SBML annotation on ListOfMembers)

    """
    def __init__(self):
        self.annotation = {}

class Group(Fbase):
    """
    Container for SBML groups

    """
    members = None
    member_ids = None
    kind = 'collection'
    _kinds_ = ('collection', 'partonomy', 'classification')
    _group_member_ids_ = None
    _member_attributes_ = None

    def __init__(self, pid):
        #delattr(self, 'compartment') # doesn't work
        pid = str(pid)
        self.setId(pid)
        self.members = []
        self.member_ids = []
        self._group_member_ids_ = []
        self.annotation = {}
        self._member_attributes_ = GroupMemberAttributes()

    def clone(self):
        """
        Return a clone of this object. Note the for Groups this is a shallow copy, in that the reference
        objects themselves are not cloned only the group (and attributes)

        """
        if type(self.__TRASH__) == dict:
            self.__TRASH__.clear()
        else:
            self.__TRASH__ = None
        cpy = copy.deepcopy(self)
        cpy.members = []
        cpy.members.extend(self.members[:])
        return cpy

    def addMember(self, obj):
        """
        Add member CBMPy object(s) to the group

         - *obj* either a single, tuple or list of CBMPy objects

        """
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        for o_ in obj:
            if isinstance(o_, Fbase):
                if o_.getId() in self.member_ids:
                    print('ERROR object {} already exist in group.\n'.format(o_.getId()))
                else:
                    self.members.append(o_)
                    self.member_ids.append(o_.getId())
                    if isinstance(o_, Group):
                        self._group_member_ids_.append(o_.getId())
            else:
                print('ERROR: object {} is not a valid CBMPy object .\n'.format(o_))

    def deleteMember(self, oid):
        """
        Deletes a group member with group id.

         - *oid* group member id

        """
        if oid in self.member_ids:
            oidx = self.member_ids.index(oid)
            obj = self.members.pop(oidx)
            self.member_ids.pop(oidx)
            if isinstance(obj, Group):
                self._group_member_ids_.remove(oid)
        else:
            print('\nError object {} is not a member of group.\n'.format(oid))

    def hasMember(self, mid):
        """
        Returns a boolean indicating whether a member is in the group.

         - *mid* the id to check

        """
        return mid in self.member_ids

    def getMembers(self, as_set=False):
        """
        Return the member objects of the group.

         - *as_set* return objects as a set rather than a list

        """
        if not as_set:
            return self.members
        else:
            return set(self.members)

    def getMemberIDs(self, as_set=False):
        """
        Return the ids of the member objects.

         - *as_set* return id's as a set rather than a list

        """
        if not as_set:
            return self.member_ids
        else:
            return set(self.member_ids)

    def getKind(self):
        """
        Return the group kind

        """
        return self.kind

    def setKind(self, kind):
        """
        Sets the kind or type of the group, this must be one of: 'collection', 'partonomy', 'classification'.

         - *kind* the kind

        """
        assert kind in self._kinds_, '\nKind must be one of: {}\n'.format(self._kinds_)
        self.kind = kind

    def addSharedMIRIAMannotation(self, qual, entity, mid):
        """
        Add a qualified MIRIAM annotation or entity to the list of members (all) rather than the group itself:

         - *qual* a Biomodels biological qualifier e.g. "is" "isEncodedBy"
         - *entity* a MIRIAM resource entity e.g. "ChEBI"
         - *mid* the entity id e.g. CHEBI:17158 or fully qualifies url (if only_qual_uri)

        """
        self._member_attributes_.addMIRIAMannotation(qual, entity, mid)

    def setSharedAnnotation(self, key, value):
        """
        Sets the list of members (all) annotation as a key : value pair.

         - *key* the annotation key
         - *value* the annotation value

        """
        self._member_attributes_.setAnnotation(key, value)

    def setSharedSBOterm(self, sbo):
        """
        Set the SBO term for the the members of the group (all).

         - *sbo* the SBOterm with format: "SBO:<7 digit integer>"

        """
        self._member_attributes_.setSBOterm(sbo)

    def setSharedNotes(self, notes):
        """
        Sets the group of objects notes attribute (all):

         - *notes* the note string, should preferably be (X)HTML for SBML

        """
        self._member_attributes_.setNotes(notes)

    def getSharedMIRIAMannotations(self):
        """
        Return a dictionary of the shared member MIRIAM annotations (rather than the group attribute).

        """
        return self._member_attributes_.getMIRIAMannotations()

    def getSharedAnnotations(self):
        """
        Return a dictionary of the shared member annotations (rather than the group attribute).

        """
        return self._member_attributes_.getAnnotations()

    def getSharedSBOterm(self):
        """
        Return the shared member SBO term (rather than the group attribute).

        """
        return self._member_attributes_.getSBOterm()

    def getSharedNotes(self):
        """
        Return the shared member notes (rather than the group attribute).

        """
        return self._member_attributes_.getNotes()

    def assignSharedAnnotationToMembers(self):
        """
        This function merges or updates the group member objects annotations with the group shared annotation.

        """
        print('INFO: Assigning shared CBMPy annotation to members, this cannot be undone.')
        if len(self._member_attributes_.annotation) > 0:
            for m_ in self.members:
                m_.annotation.update(self._member_attributes_.annotation)

    def assignSharedMIRIAMannotationToMembers(self):
        """
        This function merges or updates the group member objects MIRIAM annotations with the group shared MIRIAM annotation.

        """
        print('INFO: Assigning shared MIRIAM annotation to members, this cannot be undone.')
        if self._member_attributes_.miriam != None:
            annot = self._member_attributes_.miriam.getAllMIRIAMUris()
            for m_ in self.members:
                for k_ in annot:
                    if len(annot[k_]) > 0:
                        for u_ in annot[k_]:
                            m_.addMIRIAMuri(k_, u_)

    def assignSharedSBOtermsToMembers(self, overwrite=False):
        """
        Assigns the group shared member SBO term to the group members.

         - *overwrite* [default=False] overwrite the target SBO term if it is defined

        """
        print('INFO: Assigning shared SBOterm to members, this cannot be undone.')
        sbo = self._member_attributes_.getSBOterm()
        if sbo != None:
            for m_ in self.members:
                if m_.__sbo_term__ == None or overwrite:
                    m_.setSBOterm(sbo)

    def assignSharedNotesToMembers(self, overwrite=False):
        """
        Assigns the group shared notes to the group members.

         - *overwrite* [default=False] overwrite the target notes if they are defined

        """
        print('INFO: Assigning shared Notes to members, this cannot be undone.')
        if self._member_attributes_.notes != '':
            for m_ in self.members:
                if m_.notes == '' or overwrite:
                    m_.notes = self._member_attributes_.notes

    def assignAllSharedPropertiesToMembers(self, overwrite=False):
        """
        Assigns all group shared properties (notes, annotations, MIRIAM annotations, SBO) to the group members.

         - *overwrite* [default=False] overwrite the target notes if they are defined

        """
        self.assignSharedSBOtermsToMembers(overwrite)
        self.assignSharedNotesToMembers(overwrite)
        self.assignSharedAnnotationToMembers()
        self.assignSharedMIRIAMannotationToMembers()


class FluxBound(Fbase):
    """A reaction fluxbound"""
    reaction = None
    operation = None
    value = None
    is_bound = None
    __param__ = None

    def __init__(self, pid, reaction, operation, value):
        pid = str(pid)
        self.setId(pid)

        self.reaction = str(reaction)
        self.setValue(value)
        assert operation in ['greater', 'greaterEqual', 'less', 'lessEqual', '>=', '<=', '=', 'equal', 'E', 'G', 'L', 'GE', 'LE']
        self.operation = operation
        if self.operation in ['greater', 'greaterEqual', '>=', 'G', 'GE']:
            self.operation = 'greaterEqual'
            self.is_bound = 'lower'
        elif self.operation in ['less', 'lessEqual', '<=', 'L', 'LE']:
            self.is_bound = 'upper'
            self.operation = 'lessEqual'
        else:
            self.is_bound = 'equality'
            self.operation = 'equal'
        self.annotation = {}
        self.compartment = None
        self.__delattr__('compartment')

    def getType(self):
        """
        Returns the *type* of FluxBound: 'lower', 'upper', 'equality' or None

        """
        if self.operation in ['greater', 'greaterEqual', '>=']:
            self.is_bound = 'lower'
        elif self.operation in ['less', 'lessEqual', '<=']:
            self.is_bound = 'upper'
        elif self.operation in ['equal', 'equality', '=']:
            self.is_bound = 'equality'
        else:
            return None
        return self.is_bound

    def getReactionId(self):
        return self.reaction

    def setReactionId(self, react):
        """
        Sets the reaction attribute of the FluxBound

        """
        self.reaction = react

    def getValue(self):
        """
        Returns the current value of the attribute (input/solution)
        """
        return self.value

    def setValue(self, value):
        """
        Sets the attribute ''value''
        """
        if numpy.isreal(value):
            self.value = value
        elif numpy.isinf(value):
            self.value = value
        else:
            self.value = float(value)


class FluxBoundBase(Fbase):
    """A refactored and streamlined FluxBound base class that can be a generic bound, superclass to FluxBoundUpper and FluxBoundLower"""

    _parent = None
    operator = None
    value = None
    __param__ = None

    def __init__(self, pid, operator, value, parent=None):
        """
        - *pid* object id
        - *operator* <> GE or LE
        - *value* a float
        - *parent* [default=None] the parent reaction object

        """
        pid = str(pid)
        self.setId(pid)
        if parent is Reaction or parent is None:
            self._parent = parent
        else:
            raise RuntimeError("Invalid parent object: " + str(parent))

        if self.operator in ['greater', 'greaterEqual', '>=', 'G', 'GE']:
            self.operator = '>='
        elif self.operator in ['less', 'lessEqual', '<=', 'L', 'LE']:
            self.operator = '<='
        else:
            raise RuntimeError('Invalid operator: ' + operator)

        self.setValue(value)
        self.annotation = {}
        self.compartment = None
        #self.__delattr__('compartment')

    def getType(self):
        """
        Returns the *type* of FluxBound: 'lower', 'upper'

        """
        if self.operator  == '>=':
            return 'lower'
        else:
            return 'upper'

    def getReactionId(self):
        if self._parent is not None:
            return self._parent.getId()
        else:
            return None

    def getValue(self):
        """
        Returns the current value of the attribute (input/solution)
        """
        return self.value

    def setValue(self, value):
        """
        Sets the value attribute:

        - *value* a float

        """
        if numpy.isreal(value) or numpy.isinf(value):
            self.value = value
        else:
            print('Invalid value: ' + value)
            return False
        return True

class FluxBoundUpper(FluxBoundBase):
    def __init__(self, reaction, value=numpy.inf):
        """
        Upper Bound class, less flexible than generic superclass (no input checking) for model instantiation.

        - *value* [default=inf] a float
        - *reaction* the parent Reaction

        """
        self.setId('{}_upper_bnd'.format(reaction.getId()))
        self.operator = '<='
        self.setValue(value)
        self._parent = reaction
        self.annotation = {}
        self.compartment = None
        #self.__delattr__('compartment')


class FluxBoundLower(FluxBoundBase):
    def __init__(self, reaction, value=numpy.NINF):
        """
        Lower Bound Class, less flexible than generic superclass (no input checking) for model instantiation.

        - *value* [default=-inf] a float
        - *reaction* the parent Reaction

        """
        self.setId('{}_lower_bnd'.format(reaction.getId()))
        self.operator = '>='
        self.setValue(value)
        self._parent = reaction
        self.annotation = {}
        self.compartment = None
        #self.__delattr__('compartment')


class Parameter(Fbase):
    """Holds parameter information"""

    _associations_ = None
    constant = True
    value = None
    _is_fluxbound_ = False

    def __init__(self, pid, value, name=None, constant=True):
        """
        Parameter definition class

         - *pid* the unique parameter pid
         - *value* the value
         - *name* [default=''] the parameter name
         - *constant* [default=True] is the paramter constant (an SBML thing)

        """
        pid = str(pid)
        self.setId(pid)

        self.name = name
        self.value = value
        self.constant = constant
        self._associations_ = []
        self.annotation = {}

    def getValue(self):
        """
        Returns the current value of the attribute (input/solution)

        """
        return self.value

    def setValue(self, value):
        """
        Sets the attribute ''value''

        """
        if numpy.isreal(value):
            self.value = value
        elif numpy.isinf(value):
            self.value = value
        else:
            self.value = float(value)

    def getAssociations(self):
        """
        Return the Object ID's associated with this parameter

        """
        return self._associations_

    def addAssociation(self, assoc):
        """
        Add an object ID to associate with this object

        """
        self._associations_.append(assoc)

    def deleteAssociation(self, assoc):
        """
        Delete the object id associated with this object

        """
        if assoc in self._associations_:
            self._associations_.remove(assoc)


class Reaction(Fbase):
    """Holds reaction information"""
    reagents = None
    reversible = None
    is_exchange = False
    value = None
    reduced_cost = None
    is_balanced = None
    fva_min = None
    fva_max = None
    fva_status = None
    __bound_history__ = None
    __is_active__ = True
    _modifiers_ = None
    """
    # TODO: next major revision 0.8 is to get rid of fluxbound array
    # by adding fluxbound objects directly to the reactions this should simplify the
    # data structure but will mean a major rewrite of existing code and potentially breaking
    # backwards compatability
    upper_bound = numpy.Inf
    lower_bound = -numpy.Inf
    """
    def __init__(self, pid, name=None, reversible=True):
        pid = str(pid)
        if not self.__checkId__(pid):
            raise RuntimeError('ERROR: Invalid Id not set, \"{}\" is an invalid identifier.'.format(pid))
        self.id = pid
        self.name = name
        self.reagents = []
        self.reversible = reversible
        self.annotation = {}
        self.__TRASH__ = {}
        self.__bound_history__ = []
        self._modifiers_ = [] # reaction modifiers from SBML, read/write only

    def addReagent(self, reag):
        """
        Adds an instantiated Reagent object to the reaction

        """
        if self.__objref__ is not None:
            if reag.getId() in self.__objref__().__global_id__:
                raise RuntimeError('Duplicate obj ID detected: {}'.format(reag.getId()))
            else:
                self.__objref__().__pushGlobalId__(reag.getId(), reag)
            reag.__objref__ = self.__objref__
        self.reagents.append(reag)

    def createReagent(self, metabolite, coefficient):
        """
        Create a new reagent and add it to the reaction:

         - **metabolite** the metabolite name
         - **coefficient** the

         -- negative coefficient is a substrate
         -- positive coefficient is a product

        Will fail if a species reference already exists

        """
        assert metabolite not in self.getSpeciesIds(), '\nA reagent already refers to metabolite: %s' % metabolite
        rr = Reagent('%s_%s' % (self.getId(), metabolite), metabolite, coefficient)
        self.addReagent(rr)

    def getReagentObjIds(self):
        """
        Returns a list of the reagent id's. For the name of the reagents/metabolites use *<reaction>.getSpeciesIds()*

        """
        return [r.getId() for r in self.reagents]

    def getReagentRefs(self):
        """
        Returns a list of the reagents/metabolites

        """
        raise DeprecationWarning('DEPRECATED: please use <reaction>.getSpeciesIds')

    def getSpeciesIds(self):
        """
        Returns a list of the reagents/metabolites

        """
        return [r.species_ref for r in self.reagents]

    def getSpeciesObj(self):
        """
        Returns a list of the species objects that are reagents

        """
        if self.__objref__ == None:
            print('INFO: Only works when part of a model.')
            return []
        else:
            return [self.__objref__().getSpecies(r.species_ref) for r in self.reagents]


    def getFVAdata(self, roundnum=None, silent=True):
        """
        Returns the data generated by CBSolver.FluxVariabilityAnalysis() for this reaction as a tuple of
        (Flux, FVAmin, FVAmax, span) where span is abs(FVAmax - FVAmin). FVAmin or FVAmax is None this indicates no solution
        to that particular optimization (infeasible).

         - *roundnum* [default=None] the integer number of roundoff decimals the default is no rounding
         - *silent* [default=True] supress output to stdout

        """
        try:
            span = abs(self.fva_max - self.fva_min)
        except Exception:
            span = numpy.NaN
        out = [self.value, self.fva_min, self.fva_max, span]
        if roundnum != None:
            for f_ in range(len(out)):
                try:
                    out[f_] = round(out[f_], roundnum)
                except Exception:
                    pass

        if not silent:
            print('{}'.format(self.getId()))
            print('Flux:   {}\nFVAmin: {}\nFVAmax: {}\nSpan:   {}\n'.format(out[0], out[1], out[2], out[3]))
        return tuple(out)

    def getReagent(self, rid):
        """
        Return the one or more reagent objects which have *rid*:

         - *rid* a reagent *rid*

        """

        rgnt = [r for r in self.reagents if r.getId() == rid]
        if len(rgnt) == 0:
            return None
        elif len(rgnt) == 1:
            return rgnt[0]
        else:
            print('\nWarning multiple reagents defined for species: {}'.format(rid))
            return rgnt

    def changeId(self, pid):
        """
        Changes the Id of the reaction and updates associated FluxBounds

        """
        self.setId(pid)

    def setId(self, fid):
        """
        Sets the object Id

         - *fid* a valid c variable style id string

         Reimplements @FBase.setId()

        """

        fid = str(fid)
        if fid == self.id:
            return

        if not self.__checkId__(fid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(fid))

        oldId = self.getId()
        if self.__objref__ is not None:
            if fid not in self.__objref__().__global_id__:
                self.id = fid
                self.__objref__().__changeGlobalId__(oldId, self.id, self)
                for fb in self.__objref__().getFluxBoundsByReactionID(oldId):
                    if fb is not None:
                        fb.setReactionId(fid)
                for gpr_ in self.__objref__().gpr:
                    if gpr_.getProtein() == oldId:
                        gpr_.setProtein(fid)
                for obj in self.__objref__().objectives:
                    for fo in obj.fluxObjectives:
                        if fo.getReactionId() == oldId:
                            fo.setReactionId(fid)
            else:
                print('ERROR: setId() - object with id \"{}\" already exists ... ID *not* set.'.format(fid))
        else:
            self.id = fid
        # no matter if the reaction has a model ref or not setId resets reagent refs to new id
        for reag in self.reagents:
            reag.setId('{}_{}'.format(fid, reag.getSpecies()))

    def getValue(self):
        """
        Returns the current value of the flux.

        """
        return self.value

    def setValue(self, value):
        """
        Sets the attribute *value* in this case the flux.

        """
        self.value = value

    def getReagentWithSpeciesRef(self, sid):
        """
        Return the reagent object which refers to the *species* id. If there are multiple reagents that
        refer to the same species a list is returned.

         - *sid* the species/metabolite id

        """
        rgnt = [r for r in self.reagents if r.species_ref == sid]
        if len(rgnt) == 0:
            return None
        elif len(rgnt) == 1:
            return rgnt[0]
        else:
            print('WARNING: multiple reagents defined for species: {}'.format(sid))
            return rgnt

    def setStoichCoefficient(self, sid, value):
        """
        Sets the stoichiometric coefficient of a reagent that refers to a metabolite. Note *negative coefficients* are *substrates*
        while *positive* ones are *products*. At this point zero coefficients are not allowed

        - *sid* the species/metabolite id
        - *value* a floating point value != 0

        """
        S = self.getReagentWithSpeciesRef(sid)
        if S != None and not type(S) == list:
            S.setCoefficient(value)
        elif type(S) == list:
            raise RuntimeWarning('setStoichCoefficient({}) warning, species {} is referenced by multiple reagents: {}'.format(self.getId(), sid, [a.getId() for a in S]))
        else:
            print('ERROR: setStoichCoefficient: species {} does not exist'.format(sid))

    def getStoichiometry(self, use_names=False, altout=False):
        """
        Returns a list of (coefficient, species) pairs for this reaction

        - *use_names* [default = False] use species names rather than id's
        - *altout* [default = False] returns a dictionary [DEPRECATED]

        """
        if not use_names:
            out = [(r.getCoefficient(), r.species_ref) for r in self.reagents]
            if not altout:
                return out
            else:
                raise RuntimeError('getStoichiometry(altout=True) has been deprecated')
        else:
            out = [(r.getCoefficient(), self.__objref__().getSpecies(r.species_ref).getName()) for r in self.reagents]
            if not altout:
                return out
            else:
                raise RuntimeError('getStoichiometry(altout=True) has been deprecated')

    def getSubstrateIds(self, use_names=False):
        """
        Returns a list of the reaction substrates, species identifiers

        - *use_names* [defualt = False] use species names rather than id's

        """
        if not use_names:
            return [r.species_ref for r in self.reagents if r.getCoefficient() < 0.0]
        else:
            return [self.__objref__().getSpecies(r.species_ref).getName() for r in self.reagents if r.getCoefficient() < 0.0]

    def getProductIds(self, use_names=False):
        """
        Returns a list of the reaction products, species identifiers

        - *use_names* [default = False] use species names rather than id's

        """
        if not use_names:
            return [r.species_ref for r in self.reagents if r.getCoefficient() > 0.0]
        else:
            return [self.__objref__().getSpecies(r.species_ref).getName() for r in self.reagents if r.getCoefficient() > 0.0]

    def getGPRassociationString(self, use_labels=True):
        """
        Return the GPR string associated with this reaction (assuming it exists) or None.

        - *use_labels* [default=True] return string with lab

        """
        out = None
        try:
            out = self.__objref__().getGPRforReaction(self.getId())
            if out is not None:
                out = out.getAssociationStr(use_labels=use_labels)
        except Exception as ex:
            print(ex)
            out = None
        return out

    def deleteReagentWithSpeciesRef(self, sid):
        """
        Delete a reagent (or reagents) that refers to the species id:

         - *sid* a species/metabolite id

        """
        reags = self.getSpeciesIds()
        assert sid in reags, '\nThats not a good metabolite/species ref'
        for rr in range(len(self.reagents)-1,-1,-1):
            if self.reagents[rr].getSpecies() == sid:
                rg = self.reagents.pop(rr)
                print('Deleting reagent: {}'.format(rg.getId()))
                #removed until I have a more secure way of doing this
                #self.__TRASH__.update({rg.getId() : rg.clone()})
                if self.__objref__ is not None:
                    self.__objref__().__popGlobalId__(rg.getId())
        del rg

    #removed until I have a more secure way of doing this
    #def undeleteReagentWithSpeciesRef(self, sid):
        #"""
        #Attempts to unDelete reagent deleted with deleteReagent() that refers to the species id:

         #- *sid* a species/metabolite id

        #"""
        ##assert self.getId()+sid in self.__TRASH__, '\nAha, yes, sure, maybe, perhaps ...'
        #for rg in list(self.__TRASH__.keys()):
            #if self.__TRASH__[rg].species_ref == sid:
                #self.addReagent(self.__TRASH__.pop(rg))

    def getLowerBound(self):
        """
        Get the value of the reactions lower bound

        """
        try:
            return self.__objref__().getReactionLowerBound(self.id)
        except AttributeError as why:
            print('WARNING: This function requires that this reaction object be added to a CBMPy instance to work.')
            return None

    def getUpperBound(self):
        """
        Get the value of the reactions upper bound

        """
        try:
            return self.__objref__().getReactionUpperBound(self.id)
        except AttributeError as why:
            print('WARNING: This function requires that this reaction object be added to a CBMPy instance to work.')
            return None

    def setLowerBound(self, value):
        """
        Set the value of the reactions lower bound

         - *value* a floating point value

        """
        try:
            self.__objref__().setReactionLowerBound(self.id, value)
        except AttributeError as why:
            print('WARNING: This function requires that this reaction object be added to a CBMPy instance to work.')

    def setUpperBound(self, value):
        """
        Set the value of the reactions upper bound

         - *value* a floating point value

        """
        try:
            self.__objref__().setReactionUpperBound(self.id, value)
        except AttributeError as why:
            print('WARNING: This function requires that this reaction object be added to a CBMPy instance to work.')

    def deactivateReaction(self, lower=0.0, upper=0.0, silent=True):
        """
        Deactivates a reaction by setting its bounds to lower and upper. Restore with reactivateReaction()

         - *lower* [default=0.0] bound
         - *upper* [default=0.0] bound

        """
        self.__bound_history__ = None
        lb = self.getLowerBound()
        ub = self.getUpperBound()
        self.__bound_history__ = (lb, ub)
        self.setLowerBound(lower)
        self.setUpperBound(upper)
        self.__is_active__ = False
        if not silent:
            print('Reaction {} bounds set to [{} : {}]'.format(self.id, lower, upper))

    def reactivateReaction(self, silent=True):
        """
        Activates a reaction deactivated with deactivateReaction

        """
        if self.__bound_history__ != None:
            self.setLowerBound(self.__bound_history__[0])
            self.setUpperBound(self.__bound_history__[1])
            self.__bound_history__ = None
            self.__is_active__ = True
            if not silent:
                print('Reaction {} bounds set to [{} : {}]'.format(self.id, self.__bound_history__[0], self.__bound_history__[1]))

    def getEquation(self, reverse_symb='=', irreverse_symb='>', use_names=False):
        """
        Return a pretty printed string containing the reaction equation

         - *reverse_symb* [default = '='] the symbol to use for reversible reactions
         - *irreverse_symb* [default = '>'] the symbol to use for irreversible reactions
         - *use_names* [defualt = False] use species names rather than id's

        """
        sub = ''
        prod = ''
        for r in self.reagents:
            coeff = r.getCoefficient()
            if coeff < 0.0:
                if abs(coeff) == 1.0:
                    if not use_names:
                        sub += '{} + '.format(r.species_ref)
                    else:
                        sub += '{} + '.format(self.__objref__().getSpecies(r.species_ref).getName())
                else:
                    if not use_names:
                        sub += '({}) {} + '.format(coeff, r.species_ref)
                    else:
                        sub += '({}) {} + '.format(coeff, self.__objref__().getSpecies(r.species_ref).getName())
            else:
                if abs(coeff) == 1.0:
                    if not use_names:
                        prod += '{} + '.format(r.species_ref)
                    else:
                        prod += '{} + '.format(self.__objref__().getSpecies(r.species_ref).getName())
                else:
                    if not use_names:
                        prod += '({}) {} + '.format(coeff, r.species_ref)
                    else:
                        prod += '({}) {} + '.format(coeff, self.__objref__().getSpecies(r.species_ref).getName())
        #print(sub)
        #print(prod)
        if self.reversible:
            eq = '{} {} {}'.format(sub[:-3], reverse_symb, prod[:-2])
        else:
            eq = '{} {} {}'.format(sub[:-3], irreverse_symb, prod[:-2])
        return eq


class Species(Fbase):
    """
    Holds species/metabolite information

    """
    chemFormula = None
    charge = None
    value = None
    is_boundary = False
    reagent_of = None
    shadow_price = None

    def __init__(self, pid, boundary=False, name=None, value=float('nan'),\
                    compartment=None, charge=None, chemFormula=None):
        """
        Species/metabolite definition class

         - **pid** the unique species pid
         - **boundary** [default=False] whether the species is a variable (False) or is a boundary parameter (fixed)
         - **name** [default=''] the species name
         - **value** [default=nan] the value *not currently used*
         - **compartment** [default=None] the compartment the species is located in
         - **charge** [default=None] the species charge
         - **chemFormula** [default=None] the chemical formula

        """
        pid = str(pid)
        if not self.__checkId__(pid):
            raise RuntimeError('ERROR: Invalid Id not set, \"{}\" is an invalid identifier.'.format(pid))
        self.id = pid
        self.name = name
        self.value = value
        self.is_boundary = boundary
        self.compartment = compartment
        self.charge = charge
        self.chemFormula = chemFormula
        self.reagent_of = []
        self.annotation = {}

    def setId(self, fid, allow_rename=False):
        """
        Sets the object Id

         - *fid* a valid c variable style id string

         Reimplements @FBase.setId()

        """

        fid = str(fid)
        if fid == self.id:
            return

        if not self.__checkId__(fid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(fid))

        NEWID = False
        oldId = self.getId()
        if fid == oldId:
            return
        if self.__objref__ is not None:
            if fid not in self.__objref__().__global_id__:
                NEWID = True

            if NEWID or allow_rename:
                rids = [a[0] for a in self.__objref__().getFluxesAssociatedWithSpecies(oldId)]
                if NEWID:
                    self.id = fid
                    self.__objref__().__changeGlobalId__(oldId, self.id, self)
                for r_ in self.__objref__().reactions:
                    if r_.getId() in rids:
                        rr = r_.getReagentWithSpeciesRef(oldId)
                        if type(rr) != list:
                            rr.setSpecies(fid)
                        else:
                            for rr2 in rr:
                                rr2.setSpecies(fid)
            else:
                print('ERROR: setId() - object with id \"{}\" already exists ... ID *not* set.'.format(fid))
        else:
            self.id = fid

        if allow_rename and not NEWID:
            self.__objref__().species.pop(self.__objref__().species.index(self))

    def getValue(self):
        """
        Returns the current value of the attribute (input/solution)

        """
        return self.value

    def setValue(self, value):
        """
        Sets the attribute ''value''

        """
        self.value = value

    def isReagentOf(self):
        """
        Returns a dynamically generated list of reactions that this species occurs as a reagent

        """
        assert self.__objref__ != None, "\nWARNING: needs to be added to a model (cmod.addSpecies()) to work"
        self.reagent_of = [r.getId() for r in self.__objref__().reactions if self.getId() in r.getSpeciesIds()]
        return self.reagent_of

    def getReagentOf(self):
        """
        Returns a list of reaction id's that this metabolite occurs in

        """
        print('INFO: The static .getReagentOf() method is deprecated, please update your code to use: \".isReagentOf()\"')
        return self.isReagentOf()

    def setReagentOf(self, rid):
        """
        Adds the supplied reaction id to the reagent_of list (if it isn't one already)

         - *rid* a valid reaction id

        """
        raise RuntimeError('\nINFO: The static .setReagentOf() method is deprecated, please update your code to use: \".isReagentOf()\"')

    def setChemFormula(self, cf):
        """
        Sets the species chemical formula

         - *cf* a chemical formula e.g. CH3NO2

        """
        if cf != '' and not checkChemFormula(cf, quiet=True):
            print('ERROR: invalid chemFormula \"{}\" being set on Species \"{}\"'.format(cf, self.getId()))
        self.chemFormula = cf

    def getChemFormula(self):
        """
        Returns the species chemical formula

        """
        return self.chemFormula

    def setCharge(self, charge):
        """
        Sets the species charge:

         - *charge* a signed double but generally a signed int is used

        """
        charge = float(charge)
        self.charge = charge

    def getCharge(self):
        """
        Returns the species charge

        """
        return self.charge

    def setBoundary(self):
        """
        Sets the species so it is a boundary metabolite or fixed which does not occur in the stoichiometric matrix N

        """
        self.is_boundary = True

    def unsetBoundary(self):
        """
        Unsets the species boundary attribute so that the metabolite is free and therefore occurs in the stoichiometric matrix N

        """
        self.is_boundary = False

    ## this functionality is part of setId()
    #def rename(self, newid, overwrite=True):
        #"""
        #Changes the species id and updates all reagents in the model reactions. Note that existing species with id == newid
        #will be overwritten/deleted.

         #- *newid* the new species id.
         #- *overwrite* [default=True] overwrite species objects (highly recommended)

        #"""
        #assert self.__objref__ is not None, "\nWARNING: needs to be part of a model (cmod.addSpecies()) to work"
        #mod = self.__objref__()
        #if overwrite and newid in mod.getSpeciesIds():
            #print('INFO: overwriting existing species: {}'.format(newid))
            #mod.deleteSpecies(newid)
        #for rr in mod.getFluxesAssociatedWithSpecies(self.getId()):
                #mod.getReaction(rr[0]).getReagentWithSpeciesRef(self.getId()).setSpecies(newid)
        #self.setId(newid)


class Reagent(Fbase):
    """
    Has a reactive species id and stoichiometric coefficient:
     - negative = substrate
     - positive = product
     - species_ref a reference to a species obj

    """
    coefficient = None
    role = None
    species_ref = None
    _value_is_ref_ = False

    def __init__(self, pid, species_ref, coef):
        """
        Instantiates a reagent from a metatabolite and coefficient, note that now the coefficient
        can be a Parameter object in which case a connection is made to the linked Parameter

         - *pid* a unique id
         - *species_ref* a reference to a species id
         - *coefficient* the stoichiometric coefficient, a non-zero integer or Parameter

        """
        self.setId(pid)
        self.species_ref = species_ref
        self.setCoefficient(coef)
        self.annotation = {}
        self.compartment = None
        self.__delattr__('compartment')

    def setCoefficient(self, coef):
        """
        Sets the reagent coefficient and role, negative coefficients are substrates and positive ones are products

         - *coeff* the new coefficient

        """
        if type(coef) is Parameter:
            self.coefficient = weakref.ref(coef)
            self._value_is_ref_ = True
            value = coef.getValue()
        else:
            self.coefficient = coef
            self._value_is_ref_ = False
            value = coef

        if value < 0.0:
            self.role = 'substrate'
        elif value > 0.0:
            self.role = 'product'
        else:
            self.role = None
            print('WARNING - setCoefficient(): Zero coefficient detected {}!'.format(self.getId()))
            # raise RuntimeError('Zero coefficient detected and are currently not supported: ({}) {}!' % (value, self.getId()))


    def getCoefficient(self):
        """
        Returns the reagent coefficient

        """

        if not self._value_is_ref_:
            value = self.coefficient
        else:
            value = self.coefficient().getValue()

        return value

    def setSpecies(self, spe):
        """
        Sets the metabolite/species that the reagent reference refers to

        """
        self.species_ref = spe

    def getSpecies(self):
        """
        Returns the metabolite/species that the reagent reference refers to

        """
        return self.species_ref

    def getRole(self):
        """
        Returns the reagents role, "substrate", "product" or None

        """
        #if self.coefficient < 0.0:
            #self.role = 'substrate'
        #elif self.coefficient > 0.0:
            #self.role = 'product'
        #else:
            #self.role = None
        return self.role

## bgoli concept future reagent
#class Reagent():
    #_value = None
    #_weakref_ = False

    #@property
    #def coefficient(self):
        #print(self._value, self._weakref_)

    #@coefficient.setter
    #def coefficient(self, value):
        #if type(value) is Parameter:
            #self._value = weakref.ref(value)
            #self._weakref_ = True
        #else:
            #self._value = value
            #self._weakref_ = False

    #@coefficient.getter
    #def coefficient(self):
        #if self._weakref_:
            #x = self._value().value
        #else:
            #x = self._value
        #return x


class Gene(Fbase):
    """
    Contains all the information about a gene (or gene+protein construct depending on your philosophy)

    """
    active0 = False
    active = False
    label = None

    def __init__(self, pid, label=None, active=True):
        """
        A gene construct

         - *pid* the gene id
         - *label* the gene label this may or may not be a legal Sid
         - *active* is the gene is active or not (boolean)

        """
        pid = str(pid)
        if not self.__checkId__(pid):
            raise RuntimeError('ERROR: Invalid Id not set, \"{}\" is an invalid identifier.'.format(pid))
        self.id = pid

        if label == None:
            label = pid
            pid = 'g_' + pid

        self.setLabel(label)
        self.active0 = active
        self.active = active
        self.annotation = {}

    def setId(self, fid):
        """
        Sets the object Id

         - *fid* a valid c variable style id string

         Reimplements Fbase method

        """
        fid = str(fid)
        if fid == self.id:
            return

        if not self.__checkId__(fid):
            raise RuntimeError('ERROR: Id not set, \"{}\" is an invalid identifier.'.format(fid))

        if self.__objref__ is not None:
            if fid not in self.__objref__().__global_id__:
                mod = self.__objref__()
                gprs = mod.getGPRIdAssociatedWithGeneId(self.id)
                old_id = self.id
                mod.__changeGlobalId__(self.id, fid, self)
                self.id = fid
                if gprs is not None:
                    for gp in gprs:
                        GPR = mod.getGPRassociation(gp)
                        GPR.__renameGeneIdRefsInGPRTree__(GPR.tree, old_id, fid)
                        GPR.generefs = GPR.__getGeneRefsfromGPRDict__(GPR.tree, [])
            else:
                print('ERROR: setId() - object with id \"{}\" already exists ... ID *not* set.'.format(fid))
        else:
            self.id = fid

    def getLabel(self):
        """
        Returns the gene label

        """
        return self.label

    def setLabel(self, label):
        """
        Sets the gene label, checks that the new label is unique

        """
        if self.__objref__ is not None:
            if label in [g.getLabel() for g in self.__objref__().genes]:
                print('WARNING: setLabel() - gene with label \"{}\" already exists ... label *not* set.'.format(label))
            else:
                self.label = label
        else:
            self.label = label

    def setActive(self):
        """
        Set the gene to be active
        """
        self.active = True
        self.__objref__().__check_gene_activity__ = True

    def setInactive(self):
        """
        Set the gene to be inactive
        """
        self.active = False
        self.__objref__().__check_gene_activity__ = True

    def isActive(self):
        """
        Returns whether the gene is active or not
        """
        return self.active

    def resetActivity(self):
        """
        Reset the gene to its default activity state
        """
        self.active = self.active0
        self.__objref__().__check_gene_activity__ = True


class GeneProteinAssociation(Fbase):
    """
    This class associates genes to proteins.

    """
    #_MODIFIED_ASSOCIATION_ = False
    assoc = None
    #assoc0 = None
    protein = None
    __evalass__ = 'None'
    __evalass_compiled__ = None
    __evalass_result__ = None
    use_compiled = False
    _gene_id_ucntr_ = 0
    generefs = None
    tree = None

    def __init__(self, pid, protein, use_compiled=False):
        """
        Create a GeneProteinAssociation

         - *pid* a unique id
         - *protein* the protein the gene association referes to, in most cases this should be a reaction id
         - *use_compiled* [default=False] used compiled expressions for evaluation, potentially less portable

        """
        self.setPid(pid)
        self.generefs = []
        self.protein = protein
        self.annotation = {}
        if use_compiled:
            self.use_compiled = True

    #def evalAssociation(self):
        #"""
        #Returns an integer value representing the logical associations or None.

        #"""
        #out = None
        #_model_ = self.__objref__()
        #try:
            #out = eval(self.__evalass__)
        #except SyntaxError:
            #raise RuntimeWarning('\nError in GPR associated with reaction: %s\n%s' % (self.protein, self.assoc))
        #del _model_
        #return out

    def evalAssociation(self):
        """
        Returns an integer value representing the logical associations or None.

        """
        out = None
        if not self.use_compiled:
            try:
                out = eval(self.__evalass__)
            except SyntaxError:
                raise RuntimeWarning('\nError in GPR associated with reaction: %s\n%s' % (self.protein, self.assoc))
        else:
            try:
                eval(self.__evalass_compiled__)
                out = self.__evalass_result__
            except:
                raise RuntimeWarning('\nError in compiled GPR associated with reaction: %s\n%s' % (self.protein, self.assoc))
        return out

    def buildEvalFunc(self):
        #print(self.getTree())
        self.__evalass__ = 'int({})'.format(self.__getAssociationEvalFromGprDict__(
            self.getTree(), '', ''))
        if self.use_compiled:
            self.__evalass_compiled__ = compile('self.__evalass_result__ = {}'.format(self.__evalass__), '<inline>', 'single')

    def addGeneref(self, geneid):
        """
        Add a gene reference to the list of gene references

         - *geneid* a valid model Gene id

        """
        if geneid not in self.generefs:
            self.generefs.append(geneid)

    def deleteGeneref(self, gid):
        """
        Deletes a gene reference

        - *geneid* a valid model Gene id

        """
        if gid in self.generefs:
            self.generefs.remove(gid)

    def addAssociation(self, assoc):
        """
        Add a gene/protein association expression

        """
        #self.assoc = assoc
        raise RuntimeError('\nThis method has ceased to exist')

    def createAssociationAndGeneRefsFromTree(self, gprtree, altlabels=None):
        """
        Evaluate the GPR tree and add the genes necessary to evaluate it
        Note that this GPR should be added to a model with cmod.addGPRAssociation() before calling this method

         - *gprtree* the CBMPy GPR tree data structure
         - *altlabels* [default=None] a dictionary containing a label<-->id mapping

        """
        if self.__objref__() == None:
            raise RuntimeError("\nPlease add this GeneAssociation to a model with cmod.addGPRAssociation() before calling this method!")
        if altlabels is None:
            altlabels = {}
        genelist = self.__objref__().genes
        mod_genes = [g.getId() for g in genelist]
        react_gene = {}
        self.generefs = []
        self.setTree(gprtree)
        genes = self.__getGeneRefsfromGPRDict__(self.getTree(), [])
        for gid in genes:
            if gid in altlabels:
                label = altlabels[gid]
            else:
                label = gid
            newgid = fixId(gid, replace='_')
            if gid != newgid:
                # This needs to be tested
                newgid = fixId(gid, replace='_{}_'.format(self._gene_id_ucntr_))
                self._gene_id_ucntr_ += 1
                print('INFO: geneLabel is not Sid compatible, replacing \"{}\" with {} in geneId'.format(gid, newgid))
                #assoc = self.assoc.replace(gid, newgid)
                self.__renameGeneIdRefsInGPRTree__(self.getTree(), gid, newgid)
                self._MODIFIED_ASSOCIATION_ = True
                gid = newgid
            if gid in self.generefs:
                #print('gid in generef')
                pass
            elif gid in mod_genes:
                self.addGeneref(gid)
                #print('addGeneRef')
            else:
                #print('createAssociationAndGeneRefs', gid, label)
                if gid is None or gid is 'None':
                    #print(self.getTree(), genes, gid, label)
                    continue
                self.__objref__().addGene(Gene(gid, label, active=True))
                self.addGeneref(gid)
        self.buildEvalFunc()

    def createAssociationAndGeneRefs(self):
        raise RuntimeError("\n\nDEPRECATED CHANGE NOW!")

    def createAssociationAndGeneRefsFromString(self, assoc, altlabels=None):
        """
        Evaluate the gene/protein association and add the genes necessary to evaluate it
        Note that this GPR should be added to a model with cmod.addGPRAssociation() before calling this method

         - *assoc* the COBRA style gene protein association
         - *altlabels* [default=None] a dictionary containing a label<-->id mapping

        """
        if self.__objref__() == None:
            raise RuntimeError("\nPlease add this GeneAssociation to a model with cmod.addGPRAssociation() before calling this method!")
        genelist = self.__objref__().genes
        mod_genes = [g.getId() for g in genelist]
        if altlabels is None:
            altlabels = {}
        gene_label_id_map = {g.label:g.getId() for g in genelist}

        react_gene = {}
        self.generefs = []
        if assoc != None and assoc != '':
            #genes, self.assoc = extractGeneIdsFromString(assoc, return_clean_gpr=True)
            assoc = assoc.replace(' OR ', ' or ').replace(' AND ',' and ')
            try:
                ast.parse(assoc)
            except SyntaxError:
                print('Error in Gene Association String: {}'.format(assoc))
                old_gids, assoc2 = extractGeneIdsFromString(assoc, return_clean_gpr=True)

                old_gids.sort()
                old_gids.reverse()
                dupidhack = []
                for gid in range(len(old_gids)):
                    if old_gids[gid] in gene_label_id_map:
                        dupidhack.append(gene_label_id_map[old_gids[gid]])
                    else:
                        dupidhack.append(old_gids[gid])
                #print(old_gids)
                #print(dupidhack)

                tempids = ['{:04d}'.format(i+1) for i in range(len(old_gids))]
                rep_map = {}
                for g_ in range(len(old_gids)):
                    rep_map[tempids[g_]] = '\"{}\"'.format(dupidhack[g_])
                    assoc2 = assoc2.replace(old_gids[g_], tempids[g_])
                for id_ in tempids:
                    assoc2 = assoc2.replace(id_, rep_map[id_])
                assoc = assoc2
                del assoc2, old_gids, tempids, rep_map
            try:
                newtree = getGPRasDictFromString(ast.parse(assoc).body[0], {})
                #print(newtree)
            except SyntaxError:
                err = '\nError in Gene Association String: {}\n'.format(assoc)
                #print(err)
                raise SyntaxError(err)
            self.setTree(newtree)
            genes = self.__getGeneRefsfromGPRDict__(newtree, [])
            genes.sort()
            if len(genes) == 0:
                self.generefs = []
            else:
                for gid in genes:
                    if gid in altlabels:
                        label = altlabels[gid]
                    else:
                        label = gid
                    newgid = fixId(gid, replace='_')
                    if gid != newgid:
                        newgid = fixId(gid, replace='_{}_'.format(self._gene_id_ucntr_))
                        self._gene_id_ucntr_ += 1
                        print('INFO: geneLabel is not Sid compatible, replacing \"{}\" with {} in geneId'.format(gid, newgid))
                        #assoc = self.assoc.replace(gid, newgid)
                        self.__renameGeneIdRefsInGPRTree__(self.getTree(), gid, newgid)
                        self._MODIFIED_ASSOCIATION_ = True
                        gid = newgid
                    if gid in self.generefs:
                        #print('gid in generef')
                        pass
                    elif gid in mod_genes:
                        self.addGeneref(gid)
                        #print('addGeneRef')
                    else:
                        #print('createAssociationAndGeneRefs\n', gid, label, assoc, self.assoc)
                        self.__objref__().addGene(Gene(gid, label, active=True))
                        self.addGeneref(gid)
        else:
            self.generefs = []
        self.__objref__().__updateGeneIdx__()
        self.buildEvalFunc()

    #def buildEvalFunc(self):
        #"""
        #Builds a function which evaluates the gene expressions and evaluates to an integer using
        #the following rules:

         #- True --> 1
         #- False --> 0
         #- and --> *
         #- or --> +

        #"""
        #gids = self.getGeneIds()
        ###  print gids
        #if len(gids) > 0:
            #self.__evalass__ = self.assoc
            #_model_ = self.__objref__()
            ## this is to avoid substring replacements
            #gids = sorted(gids, key=len)
            #gids.reverse()
            #for g in gids:
                #self.__evalass__ = self.__evalass__.replace(g, "_model_.genes[{}].isActive()".format(_model_.__genes_idx__.index(g)))
            #self.__evalass__ = self.__evalass__.replace(' or ', ' + ')
            #self.__evalass__ = self.__evalass__.replace(' OR ', ' + ')
            #self.__evalass__ = self.__evalass__.replace(' and ', ' * ')
            #self.__evalass__ = self.__evalass__.replace(' AND ', ' * ')
            #self.__evalass__ = 'int(%s)' % self.__evalass__
            #del _model_
            ##self.__evalass__ = compile(self.__evalass__, 'GeneAss', 'exec')

    def getGenes(self):
        """
        Return a list of gene objects associated with this GPRass
        """
        return [self.__objref__().getGene(g) for g in self.generefs]

    def getGeneLabels(self):
        """
        Return a list of gene labels associated with this GPRass
        """
        return [self.__objref__().getGene(g).getLabel() for g in self.generefs]

    def getGene(self, gid):
        """
        Return a gene object with id

        """
        if gid in self.generefs:
            return self.__objref__().getGene(gid)
        else:
            print('WARNING: {} is not a valid gene id'.format(gid))
            return None

    def getAssociationStr(self, use_labels=False):
        """
        Return the gene association string, alternatively return string with labels

        - *use_lablels* [default=False] return the gene association string with labels rather than geneId's (FBCv2 issue)

        """
        out = self.__getAssociationStrFromGprDict__(self.getTree(), '', parent='')
        if use_labels:
            keymap = {}
            for g in self.generefs:
                keymap[g] = self.__objref__().getGene(g).getLabel()
            keys = list(keymap.keys())
            keys.sort(reverse=True)
            for k in keys:
                out = out.replace(k, keymap[k])
        return out

    def __getAssociationStrFromGprDict__(self, gprd, out, parent=''):
        """
        Get a old school GPR association string from a CBMPy gprDict, e.g. obtained from gpr.getTree()

         - *gprd* the gprDictionary
         - *out* the output string
         - *parent* [default=''] the string representing the current nodes parent relationship, used for recursion

        """
        out2 = '('
        for k in gprd:
            if k.startswith('_AND_'):
                out2 += '{}{}'.format(self.__getAssociationStrFromGprDict__(gprd[k], out, ' and '), parent)
            elif k.startswith('_OR_'):
                out2 += '{}{}'.format(self.__getAssociationStrFromGprDict__(gprd[k], out, ' or '), parent)
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

    def __getGeneRefsfromGPRDict__(self, gprd, out):
        """
        Extract the gene id references from the GPR tree

        - *gprd* the gprTree
        - *out* the output list

        """
        for k in gprd:
            if k.startswith('_AND_'):
                self.__getGeneRefsfromGPRDict__(gprd[k], out)
            elif k.startswith('_OR_'):
                self.__getGeneRefsfromGPRDict__(gprd[k], out)
            else:
                if k not in out:
                    out.append(k)
        return out

    def __renameGeneIdRefsInGPRTree__(self, gprd, old, new):
        """
        Rename gene ids in the gpr Tree, works inplace

        - *gprd* the gprTree
        - *old* the old gene id
        - *new* the new gene id

        """
        for k in list(gprd):
            if k.startswith('_AND_'):
                self.__renameGeneIdRefsInGPRTree__(gprd[k], old, new)
            elif k.startswith('_OR_'):
                self.__renameGeneIdRefsInGPRTree__(gprd[k], old, new)
            elif k == old:
                gprd.pop(k)
                gprd[new] = new

    def getGeneIds(self):
        """
        Return a list of gene id's
        """
        return self.generefs

    def getActiveGenes(self):
        """
        Return a list of active gene objects

        """
        out = []
        for g_ in self.generefs:
            G = self.__objref__().getGene(g_)
            if G.isActive():
                out.append(G)
        return out

    def getProtein(self):
        """
        Return the protein associated with this set of genes

        """
        return self.protein

    def setProtein(self, protein):
        """
        Sets the protein associated with this set of genes

        """
        self.protein = protein

    def setGeneInactive(self, gid):
        """
        Set a gene to be inactive
        """
        try:
            self.getGene(gid).setInactive()
        except Exception as ex:
            print('\nGene {} does not exist'.format(gid))

    def setGeneActive(self, gid):
        """
        Set a gene to be inactive
        """
        try:
            self.getGene(gid).setActive()
        except Exception as ex:
            print('\nGene {} does not exist'.format(gid))

    def setAllGenesActive(self):
        """
        Activate all genes in association

        """
        for g_ in self.generefs:
            G = self.__objref__().getGene(g_)
            G.setActive()

    def setAllGenesInactive(self):
        """
        Deactivates all genes in association

        """
        for g_ in self.generefs:
            G = self.__objref__().getGene(g_)
            G.setInactive()

    def isProteinActive(self):
        """
        This returns a boolean which indicates the result of evaluating the gene association. If the result is positive
        then the protein is expressed and *True* is returned, otherwise if the expression evaluates to a value of 0 then
        the protein is not expressed and  *False* is returned.
        """
        res = self.evalAssociation()
        if res == 0:
            return False
        else:
            return True

    def setTree(self, tree):
        """
        Add a GPR dictionary/tree representation to the GPR.

         - *tree* a dictionary representation of a GPR.

        """
        self.tree = tree

    def getTree(self):
        """
        Return the dictionary/tree representation of the GPR

        """
        return self.tree

    def getTreeCopy(self):
        """
        Return a copy of the dictionary/tree representation of the GPR

        """
        return copy.deepcopy(self.tree)


    def __getAssociationEvalFromGprDict__(self, gprd, out, parent=''):
        """
        Get a GPR evaluation string from a CBMPy gprDict, e.g. obtained from gpr.getTree()

         - *gprd* the gprDictionary
         - *out* the output string
         - *parent* [default=''] the string representing the current nodes parent relationship, used for recursion
         - *model* an FBA model with gene information

        """

        out2 = '('
        for k in gprd:
            if k.startswith('_AND_'):
                out2 += '{}{}'.format(self.__getAssociationEvalFromGprDict__(gprd[k], out, ' * '), parent)
            elif k.startswith('_OR_'):
                out2 += '{}{}'.format(self.__getAssociationEvalFromGprDict__(gprd[k], out, ' + '), parent)
            else:
                out2 += 'self.__objref__().getGene(\'{}\').isActive(){}'.format(k, parent)
        if out2.endswith(' + '):
            out2 = out2[:-3]
        elif out2.endswith(' * '):
            out2 = out2[:-3]
        out2 += ')'
        out = out + out2
        if ' + ' not in out and ' * ' not in out:
            out = out[1:-1]
        return out

    def deleteGeneFromAssociation(self, gid):
        """
        Deletes a gene id from the gene association. *WARNING* this process is irreversible!!

        - *gid* a valid gene identifier (not label)

        """
        #print('DeleteGene is processing gene: {}'.format(gid))
        if gid in self.generefs:
            #print('DeleteGene is deleting: {} ...'.format(gid))
            #print(self.tree)
            self.deleteGeneref(gid)
            self.tree = self.__deleteGeneFromTree__(self.tree, gid)
            self.tree = self.__pruneTree__(self.tree)
            #print(self.tree)
            self.buildEvalFunc()
        else:
            print('Gene Id {} is not part of GPR {}'.format(gid, self.getId()))

    def __deleteGeneFromTree__(self, D, delid):
        """
        Recursively delete a gene Id from a gprTree. This is a newer refactored version of the previous method
        that requires a call to __

        """
        for k in list(D):
            if k.startswith('_AND_') or k.startswith('_OR_'):
                D[k] = self.__deleteGeneFromTree__(D[k], delid)
                if len(D[k]) == 0:
                    D.pop(k)
                elif len(D[k]) == 1:
                    D.update(D.pop(k))
            elif k == delid:
                D.pop(k)
            elif len(k) == 0:
                D.pop(k)
            elif len(D[k]) == 1:
                D.update(D.pop(k))
        return D

    # old
    #def __deleteGeneFromTree__(self, D, delid):
        #"""
        #Recursively delete a gene Id from a gprTree.

        #"""
        #for k in list(D):
            #if k == delid:
                #D.pop(k)
            #elif len(k) == 0:
                #D.pop(k)
            #elif len(D[k]) == 1:
                #D.update(D.pop(k))
            #elif k.startswith('_AND_') or k.startswith('_OR_'):
                #D[k] = self.__deleteGeneFromTree__(D[k], delid)
                #if len(D[k]) == 0:
                    #D.pop(k)
                #elif len(D[k]) == 1:
                    #D.update(D.pop(k))
        #return D

    def __pruneTree__(self, D):
        """
        Recursively checks the tree for the correct number of children

        """
        for k in list(D):
            if k.startswith('_AND_') and len(D[k]) == 1:
                D.update(D.pop(k))
                print('Pruning branch: {}'.format(D))
            elif k.startswith('_OR_') and len(D[k]) == 1:
                D.update(D.pop(k))
                print('Pruning branch: {}'.format(D))
            elif k.startswith('_AND_'):
                self.__pruneTree__(D[k])
            elif k.startswith('_OR_'):
                self.__pruneTree__(D[k])
        return D


