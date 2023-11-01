import os, json, numpy, re, time
import libsbml
import libsbml as SB
#import cbmpy

#print('libSBMLversion:', libsbml.getLibSBMLVersionString())
#cmod = cbmpy.readSBML3FBC('cbmpy_test_core')

RE_loKVD1 = re.compile('<data.*?/>', re.DOTALL)
RE_loKVD1_id = re.compile('id=".*?"')
# RE_loKVD1_type = re.compile('type=".*?"')
RE_loKVD1_value = re.compile('value=".*?"')

from html import parser
HTMLParser = parser.HTMLParser



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
from html import unescape
__tagStripper__.unescape = unescape

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


print('\n##########\n\nThis is a rough FBCv2 to FBCv3 converter to test CBMPyV3 support. Only the core model components are converted ... YOU HAVE BEEN WARNED!!!\n\n##########\n')

# FBCv2init
reader = libsbml.SBMLReader()
writer = libsbml.SBMLWriter()
#D = reader.readSBML('e_coli_core.xml')
D = reader.readSBML('test_user_def_rc.xml')
with open('test_user_def_rc.xml.user_constraints.json', 'r') as F:
    UC = json.load(F)
print(UC)
assert D.getNumErrors() == 0, 'SBML file read error'
print(D)
M = D.getModel()
print(M)
Mfbc = M.getPlugin('fbc')
print(Mfbc)

# FBCv3init
fbc3ns = 'http://www.sbml.org/sbml/level3/version1/fbc/version3'
sbmlns3 = libsbml.SBMLNamespaces(3, 1, "fbc", 3)
sbmlns3.addNamespace('http://www.w3.org/1999/xhtml', 'html')
DF3 = libsbml.SBMLDocument(sbmlns3)
DF3.setPackageRequired("fbc", False)
MF3 = DF3.createModel()
MF3fbc = MF3.getPlugin('fbc')
MF3fbc.setStrict(True)


writer.writeSBMLToFile(DF3, 'the_first.fbc3.xml')


# compartments
for c_ in range(M.getNumCompartments()):
    C = M.getCompartment(c_)
    C3 = MF3.createCompartment()
    C3.setId(C.getId())
    C3.setName(C.getName())
    C3.setSize(C.getSize())
    C3.setConstant(True)

# species
for s_ in range(M.getNumSpecies()):
    S = M.getSpecies(s_)
    S3 = MF3.createSpecies()
    S3.setId(S.getId())
    S3.setBoundaryCondition(S.getBoundaryCondition())
    S3.setCompartment(S.getCompartment())
    S3.setConstant(S.getConstant())
    S3.setName(S.getName())
    S3.setInitialConcentration(S.getInitialConcentration())
    S3.setInitialAmount(S.getInitialAmount())
    S3.setHasOnlySubstanceUnits(S.getHasOnlySubstanceUnits())
    # V3 changes
    Sfbc = S.getPlugin('fbc')
    S3fbc = S3.getPlugin('fbc')
    S3fbc.setCharge(Sfbc.getChargeAsDouble())
    S3fbc.setChemicalFormula(Sfbc.getChemicalFormula())
    # spread some kv happiness around
    kvp = S3fbc.createKeyValuePair()
    kvp.setKey('time')
    kvp.setValue(time.strftime('%H:%M:%S'))




# reactions V3
for r_ in range(M.getNumReactions()):
    R =  M.getReaction(r_)
    Rfbc = R.getPlugin('fbc')

    R3 = MF3.createReaction()
    R3.setFast(False)
    R3.setReversible(True)
    R3fbc = R3.getPlugin('fbc')

    R3.setId(R.getId())
    R3fbc.setUpperFluxBound(Rfbc.getUpperFluxBound())
    R3fbc.setLowerFluxBound(Rfbc.getLowerFluxBound())

    if Rfbc.getUpperFluxBound() not in [p.getId() for p in MF3.parameters]:
        ub = MF3.createParameter()
    else:
        ub = MF3.getParameter(Rfbc.getUpperFluxBound())
    ub.setId(Rfbc.getUpperFluxBound())
    ub.setValue(M.getParameter(Rfbc.getUpperFluxBound()).getValue())
    ub.setConstant(True)

    if Rfbc.getLowerFluxBound() not in [p.getId() for p in MF3.parameters]:
        lb = MF3.createParameter()
    else:
        lb = MF3.getParameter(Rfbc.getLowerFluxBound())
    lb.setId(Rfbc.getLowerFluxBound())
    lb.setValue(M.getParameter(Rfbc.getLowerFluxBound()).getValue())
    lb.setConstant(True)

    # reactants
    for re_ in range(R.getNumReactants()):
        RE = R.getReactant(re_)
        RE3 = R3.createReactant()
        RE3.setConstant(True)
        RE3.setId(RE.getId())
        RE3.setName(RE.getName())
        RE3.setSpecies(RE.getSpecies())
        RE3.setStoichiometry(RE.getStoichiometry())

    # products
    for re_ in range(R.getNumProducts()):
        PRE = R.getProduct(re_)
        PRE3 = R3.createProduct()
        PRE3.setConstant(True)
        PRE3.setId(PRE.getId())
        PRE3.setName(PRE.getName())
        PRE3.setSpecies(PRE.getSpecies())
        PRE3.setStoichiometry(PRE.getStoichiometry())

    # KeyValuePairs
    kv_pairs = sbml_readKeyValueDataAnnotation(R.getAnnotationString())
    print(kv_pairs)

    for kv_ in kv_pairs:
        kvp = R3fbc.createKeyValuePair()
        kvp.setKey(kv_)
        if kv_pairs[kv_] is not None:
            kvp.setValue(kv_pairs[kv_])






FBC_VARIABLE_TYPES = {libsbml.FBC_VARIABLE_TYPE_LINEAR : 'linear',
                      libsbml.FBC_VARIABLE_TYPE_QUADRATIC : 'quadratic',
                      libsbml.FBC_VARIABLE_TYPE_INVALID : 'invalid'
                      }

objectives = Mfbc.getListOfObjectives()

for o_ in objectives:
    print('obj', o_.getId())
    O3 = MF3fbc.createObjective()
    O3.setId(o_.getId())
    O3.setType(o_.getType())
    for fo_ in o_.getListOfFluxObjectives():
        print(fo_)
        print(fo_.getId(), fo_.getCoefficient(), fo_.getReaction())
        fo3 = O3.createFluxObjective()
        fo3.setId(fo_.getId())
        fo3.setCoefficient(fo_.getCoefficient())
        fo3.setReaction(fo_.getReaction())
        fo3.setVariableType(libsbml.FBC_VARIABLE_TYPE_LINEAR)
        print('variable type', FBC_VARIABLE_TYPES[fo3.getVariableType()])

MF3fbc.getListOfObjectives().setActiveObjective(objectives.getActiveObjective())

# create an additional dummy quadratic objective from reactions 4 and 5
O3 = MF3fbc.createObjective()
O3.setId('DummyObjectiveToTestQuadraticCoefficients')
O3.setType('minimize')
f1 = O3.createFluxObjective()
f2 = O3.createFluxObjective()

f1.setId('dummy1')
f1.setCoefficient(4)
f1.setReaction(MF3.reactions[3].getId())
f1.setVariableType(libsbml.FBC_VARIABLE_TYPE_QUADRATIC)
f2.setId('dummy2')
f2.setCoefficient(-2)
f2.setReaction(MF3.reactions[4].getId())
f2.setVariableType(libsbml.FBC_VARIABLE_TYPE_QUADRATIC)


for u in UC:
    U = MF3fbc.createUserDefinedConstraint()
    U.setId(u)
    for f in UC[u]['fluxes']:
        print(f[0], f[1], 'linear')
        UCFO = U.createUserDefinedConstraintComponent()
        UCFO.setId('{}_{}'.format(u, f[1]))
        UCFO.setVariable(f[1])
        UCFO.setCoefficient(f[0])
        UCFO.setVariableType('linear')

    ub = MF3.createParameter()
    ub.setId('uc_{}_{}_{}'.format(u, f[1], 'ub'))
    ub.setConstant(True)

    lb = MF3.createParameter()
    lb.setId('uc_{}_{}_{}'.format(u, f[1], 'lb'))
    lb.setConstant(True)

    if UC[u]['operator'] == 'E':
        lb.setValue(UC[u]['rhs'])
        ub.setValue(UC[u]['rhs'])
    elif UC[u]['operator'] in ['G', 'GE']:
        lb.setValue(UC[u]['rhs'])
        ub.setValue(numpy.inf)
    elif UC[u]['operator'] in ['L', 'LE']:
        lb.setValue(numpy.NINF)
        ub.setValue(UC[u]['rhs'])

    U.setLowerBound(lb.getId())
    U.setUpperBound(ub.getId())






writer.writeSBMLToFile(DF3, 'the_second.fbc3.xml')
