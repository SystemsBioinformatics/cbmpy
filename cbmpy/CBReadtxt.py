"""
CBMPy: CBReadtxt module
=======================
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
Last edit: $Author: bgoli $ ($Id: CBReadtxt.py 710 2020-04-27 14:22:34Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import

# from __future__ import unicode_literals


from .CBConfig import __CBCONFIG__ as __CBCONFIG__

__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

import os
import time
import re
import numpy
from . import CBModel


# userspace
SYMB_SPLIT = ';'
COMPARTMENTS = ['[c]', '[e]', '[p]']
REPLACE_COMPARTMENTS = ['_c', '_e', '_p']
SYMB_IRR = '-->'
SYMB_REV = '<==>'
BIG_BOUND = 999999.0
# userspace


def getBounds(fname, reaction_prefix='R_', has_header=False):
    bndsF = open(fname, 'r')
    bndsLst = []
    bndsOut = {}
    for L in bndsF:
        L = [e.strip() for e in L.split(SYMB_SPLIT)]
        # general replacements '\"' --> '' and '(e)' --> '_e'
        L = [e.replace('\"', '').replace('(e)', '__e') for e in L]
        # print L
        bndsLst.append(L)
    if has_header:
        print(bndsLst.pop(0))
    for b in bndsLst:
        # note the '-' --> '_' and '.' --> '_' replacement in the id
        bndsOut.update(
            {
                reaction_prefix
                + b[0]
                .replace('-', '_')
                .replace('.', '_'): {
                    'min': float(b[1].replace(',', '.')),
                    'max': float(b[2].replace(',', '.')),
                }
            }
        )
    bndsF.close()
    return bndsOut


def getReactions(
    fname,
    reaction_prefix='R_',
    has_header=False,
    save_rpt=False,
    ignore_duplicates=False,
):
    assert os.path.exists(fname), '\n%s is not a valid file!'
    rF = open(fname, 'r')
    rxnLst = []
    Reac = {}
    dupcntr = {}
    real_duplicates = []
    gene_ors = []
    for L in rF:
        L = [e.strip() for e in L.split(SYMB_SPLIT)]
        L = [e.replace('\"', '').replace('(e)', '__e') for e in L]
        rxnLst.append(L)
    if has_header:
        print('Deleting header row:', rxnLst.pop(0))
        time.sleep(2)
    for r in rxnLst:
        assert r[0] != '', '\nReaction with no name: %s' % r
        # print r
        # note the '-' --> '_' and '.' --> '_' replacement in the id
        if len(r) >= 4:
            gene = r[3]
        else:
            gene = 'unknown'
        R_id = reaction_prefix + r[0].replace('-', '_').replace('.', '_')
        ADD_REACTION = True
        if R_id in tuple(Reac):
            if (
                r[1] == Reac[R_id]['eqn']
                and r[2] == Reac[R_id]['name']
                and gene == Reac[R_id]['gene']
            ):
                print('Complete duplicate skipping reaction', R_id)
                if R_id not in real_duplicates:
                    real_duplicates.append(R_id)
                ADD_REACTION = False
            # elif r[1] == Reac[R_id]['eqn'] and r[2] == Reac[R_id]['name']:
            elif r[1] == Reac[R_id]['eqn']:
                Reac[R_id]['gene'] = Reac[R_id]['gene'] + ' OR %s' % gene.strip()
                print(
                    'Reaction differs only by gene association adding OR association',
                    R_id,
                    Reac[R_id]['gene'],
                )
                ADD_REACTION = False
                if R_id not in gene_ors:
                    gene_ors.append(R_id)
                if r[2] not in Reac[R_id]['name']:
                    Reac[R_id]['name'] = '%s OR %s' % (Reac[R_id]['name'], r[2])
            else:
                if ignore_duplicates:
                    ADD_REACTION = False
                    print('Ignoring duplicate reaction creating new reaction:', R_id)
                else:
                    if R_id in dupcntr:
                        dupcntr[R_id] = dupcntr[R_id] + 1
                    else:
                        dupcntr.update({R_id: 1})
                    R_id = '%s_%i' % (R_id, dupcntr[R_id])
                    ADD_REACTION = True
                    print('Duplicate reaction creating new reaction:', R_id)

        if ADD_REACTION:
            ##  Reac.update({R_id : {'id' : R_id, 'name' : r[1], 'eqn' : r[2], 'gene' : gene} })
            Reac.update({R_id: {'id': R_id, 'name': r[2], 'eqn': r[1], 'gene': gene}})
    rF.close()

    if save_rpt:
        F = open('csv_read_rpt.txt', 'w')
        F.write('\nCSV Read Report\n********************\n')
        F.write('\nDuplicate entries ignored:\n')
        for d in real_duplicates:
            F.write('\t%s\n' % d)
        F.write('\nDuplicate reactions renamed (instances):\n')
        for d in dupcntr:
            F.write('\t%s : %i\n' % (d, dupcntr[d] + 1))
        if len(dupcntr) < 1:
            F.write('\tNone\n')
        F.write('\nDuplicate reactions used for gene association:\n')
        for d in gene_ors:
            F.write('\n\t%s\n\t%s\n\t%s\n' % (d, Reac[d]['gene'], Reac[d]['name']))
        F.close()

    return Reac


def getReactions_old_format(
    fname,
    reaction_prefix='R_',
    has_header=False,
    save_rpt=False,
    ignore_duplicates=False,
):
    assert os.path.exists(fname), '\n%s is not a valid file!'
    rF = open(fname, 'r')
    rxnLst = []
    Reac = {}
    dupcntr = {}
    real_duplicates = []
    gene_ors = []
    for L in rF:
        L = [e.strip() for e in L.split(SYMB_SPLIT)]
        L = [e.replace('\"', '').replace('(e)', '__e') for e in L]
        rxnLst.append(L)
    if has_header:
        print('Deleting header row:', rxnLst.pop(0))
        time.sleep(2)
    for r in rxnLst:
        assert r[0] != '', '\nReaction with no name: %s' % r
        # print r
        # note the '-' --> '_' and '.' --> '_' replacement in the id
        if len(r) >= 4:
            gene = r[3]
        else:
            gene = 'unknown'
        R_id = reaction_prefix + r[0].replace('-', '_').replace('.', '_')
        ADD_REACTION = True
        if R_id in tuple(Reac):
            if (
                r[1] == Reac[R_id]['eqn']
                and r[2] == Reac[R_id]['name']
                and gene == Reac[R_id]['gene']
            ):
                print('Complete duplicate skipping reaction', R_id)
                if R_id not in real_duplicates:
                    real_duplicates.append(R_id)
                ADD_REACTION = False
            # elif r[1] == Reac[R_id]['eqn'] and r[2] == Reac[R_id]['name']:
            elif r[1] == Reac[R_id]['eqn']:
                Reac[R_id]['gene'] = Reac[R_id]['gene'] + ' OR %s' % gene.strip()
                print(
                    'Reaction differs only by gene association adding OR association',
                    R_id,
                    Reac[R_id]['gene'],
                )
                ADD_REACTION = False
                if R_id not in gene_ors:
                    gene_ors.append(R_id)
                if r[2] not in Reac[R_id]['name']:
                    Reac[R_id]['name'] = '%s OR %s' % (Reac[R_id]['name'], r[2])
            else:
                if ignore_duplicates:
                    ADD_REACTION = False
                    print('Trying to create duplicate reaction, ignoring:', R_id)
                else:
                    if R_id in dupcntr:
                        dupcntr[R_id] = dupcntr[R_id] + 1
                    else:
                        dupcntr.update({R_id: 1})
                    R_id = '%s_%i' % (R_id, dupcntr[R_id])
                    ADD_REACTION = True
                    print('Duplicate reaction, creating new reaction:', R_id)

        if ADD_REACTION:
            ##  Reac.update({R_id : {'id' : R_id, 'name' : r[1], 'eqn' : r[2], 'gene' : gene} })
            Reac.update({R_id: {'id': R_id, 'name': r[1], 'eqn': r[2], 'gene': gene}})
    rF.close()

    if save_rpt:
        F = open('csv_read_rpt.txt', 'w')
        F.write('\nCSV Read Report\n********************\n')
        F.write('\nDuplicate entries ignored:\n')
        for d in real_duplicates:
            F.write('\t%s\n' % d)
        F.write('\nDuplicate reactions renamed (instances):\n')
        for d in dupcntr:
            F.write('\t%s : %i\n' % (d, dupcntr[d] + 1))
        if len(dupcntr) < 1:
            F.write('\tNone\n')
        F.write('\nDuplicate reactions used for gene association:\n')
        for d in gene_ors:
            F.write('\n\t%s\n\t%s\n\t%s\n' % (d, Reac[d]['gene'], Reac[d]['name']))
        F.close()

    return Reac


def parseReactions(sRxns):
    Reactions = {}
    for RXN in sRxns:
        R = sRxns[RXN]['eqn']
        Rcmp = None
        Reversible = None
        Exchange = False
        if sRxns[RXN]['id'][-3:] == '__e':
            Exchange = True
        if SYMB_REV in R:
            Reversible = True
            LHS, RHS = R.split(SYMB_REV)
        elif SYMB_IRR in R:
            Reversible = False
            LHS, RHS = R.split(SYMB_IRR)
        else:
            print('No reversability information in reaction: %s' % RXN)
            print(R)
            raw_input('\nPress any key to continue ...\n')
        if LHS[:3] in COMPARTMENTS:
            Rcmp = R[:3]
            cmp, LHS = [i.strip() for i in LHS.split(':')]

            if len(LHS) > 0:
                LHS = [i.strip().replace('-', '_') for i in LHS.split('+')]
                LHS = [i + REPLACE_COMPARTMENTS[COMPARTMENTS.index(cmp)] for i in LHS]
            if len(RHS) > 0:
                RHS = [i.strip().replace('-', '_') for i in RHS.split('+')]
                RHS = [i + REPLACE_COMPARTMENTS[COMPARTMENTS.index(cmp)] for i in RHS]
        else:
            LHS = [r.strip().replace('-', '_') for r in LHS.split('+')]
            RHS = [r.strip().replace('-', '_') for r in RHS.split('+')]
            for c in range(len(COMPARTMENTS)):
                LHS = [r.replace(COMPARTMENTS[c], REPLACE_COMPARTMENTS[c]) for r in LHS]
                RHS = [r.replace(COMPARTMENTS[c], REPLACE_COMPARTMENTS[c]) for r in RHS]

        cf = re.compile("\(.*\)")
        LHS = [[re.findall(cf, l), l] for l in LHS]
        RHS = [[re.findall(cf, r), r] for r in RHS]

        for l in range(len(LHS)):
            if LHS[l][1] == '_c':
                raw_input('Ooops1', LHS)
            if len(LHS[l][0]) == 0:
                LHS[l][0] = 1.0
                LHS[l][1] = 'M_' + LHS[l][1].replace('.', '_')
            else:
                LHS[l][1] = 'M_' + LHS[l][1].replace(LHS[l][0][0], '').strip().replace(
                    '.', '_'
                )
                LHS[l][0] = float(LHS[l][0][0].replace('(', '').replace(')', ''))
            LHS[l] = tuple(LHS[l])

        for r in range(len(RHS)):
            if RHS[r][1] == '_c':
                raw_input('Ooops2', RHS)
            if len(RHS[r][0]) == 0:
                RHS[r][0] = 1.0
                RHS[r][1] = 'M_' + RHS[r][1].replace('.', '_')
            else:
                RHS[r][1] = 'M_' + RHS[r][1].replace(RHS[r][0][0], '').strip().replace(
                    '.', '_'
                )
                RHS[r][0] = float(RHS[r][0][0].replace('(', '').replace(')', ''))
            RHS[r] = tuple(RHS[r])
        geneAss = 'unknown'
        if 'gene' in sRxns[RXN]:
            geneAss = sRxns[RXN]['gene']
        Reactions.update(
            {
                RXN: {
                    'id': RXN,
                    'reversible': Reversible,
                    'exchange': Exchange,
                    'substrates': LHS,
                    'products': RHS,
                    'name': sRxns[RXN]['name'],
                    'gene': geneAss,
                }
            }
        )
    return Reactions


def addBoundsToReactions(Reactions, Bounds, default=1111.0):
    default = abs(default)
    for r in Reactions:
        if r in Bounds:
            Reactions[r].update({'min': Bounds[r]['min'], 'max': Bounds[r]['max']})
        else:
            # print '\t:'
            if Reactions[r]['reversible']:
                Reactions[r].update({'min': -default, 'max': default})
            else:
                Reactions[r].update({'min': 0.0, 'max': default})


def dumpReactionsToTxt(Reactions, fname):
    dF = open(fname, 'w')
    for rx in Reactions:
        dF.write('\n%s\n' % Reactions[rx]['id'])
        dF.write('Name: %s\n' % Reactions[rx]['name'])
        dF.write('Reversible: %s\n' % Reactions[rx]['reversible'])
        dF.write('Exchange: %s\n' % Reactions[rx]['exchange'])
        dF.write(
            'Bounds: min=%s max=%s\n' % (Reactions[rx]['min'], Reactions[rx]['max'])
        )
        if Reactions[rx]['reversible']:
            rsign = '<==>'
        else:
            rsign = '-->'

        st = ''
        for s in Reactions[rx]['substrates']:
            st += '%s %s + ' % s
        st = st[:-2]
        pt = ''
        for p in Reactions[rx]['products']:
            pt += '%s %s + ' % p
        pt = pt[:-2]
        dF.write('%s %s %s\n' % (st.strip(), rsign.strip(), pt.strip()))
    dF.flush()
    dF.close()


def getSpecies(reactions):
    species = {}
    for r in reactions:
        for s in reactions[r]['substrates']:
            species.update({s[1]: {'compartment': s[1][-1]}})
        for p in reactions[r]['products']:
            species.update({p[1]: {'compartment': p[1][-1]}})
    return species


def dumpSpeciesToTxt(species, fname):
    dF = open(fname, 'w')
    for s in species:
        dF.write('\n%s\n' % s)
        dF.write('Compartment: %s\n' % species[s]['compartment'])
    dF.flush()
    dF.close()


def buildFBAobj(reactions, species, objfname, modname='model'):
    fm = CBModel
    M = fm.Model(modname)
    M.setName(modname)
    M.description = modname
    for s in species:
        S = fm.Species(s, compartment=species[s]['compartment'])
        S.value = 0.0
        S.setName(s)
        M.addSpecies(S)
    for r in reactions:
        R = fm.Reaction(
            r, name=reactions[r]['name'], reversible=reactions[r]['reversible']
        )
        for su in reactions[r]['substrates']:
            SU = fm.Reagent(r + su[1], su[1], -su[0])
            R.addReagent(SU)
            spx = M.getSpecies(su[1])
            if su[1] not in spx.isReagentOf():
                spx.reagent_of.append(R.getId())
        for pr in reactions[r]['products']:
            PR = fm.Reagent(r + pr[1], pr[1], pr[0])
            R.addReagent(PR)
            spx = M.getSpecies(pr[1])
            if pr[1] not in spx.isReagentOf():
                spx.reagent_of.append(R.getId())
        R.is_exchange = reactions[r]['exchange']
        LB = fm.FluxBound(r + 'min', r, 'greaterEqual', reactions[r]['min'])
        UB = fm.FluxBound(r + 'max', r, 'lessEqual', reactions[r]['max'])
        R.value = 0.0
        if 'gene' in reactions[r]:
            R.setAnnotation('GENE ASSOCIATION', reactions[r]['gene'])
        M.addReaction(R, create_default_bounds=False)
        M.addFluxBound(LB)
        M.addFluxBound(UB)
    if objfname != None:
        OB = fm.Objective('objective1', 'maximize')
        FO = fm.FluxObjective(objfname + '_obj', objfname)
        OB.addFluxObjective(FO)
        M.addObjective(OB, active=True)
    return M


def readCSV(
    model_file,
    bounds_file=None,
    biomass_flux=None,
    model_id='FBAModel',
    reaction_prefix='R_',
    has_header=False,
):
    """
    This function loads a CSV file and translates it into a Python object::

     - *model_file* the name of the CSV file that contains the model
     - *bounds_file* the name of the CSV file that contains the flux bounds
     - *biomass_flux* the name of the reaction that is the objective function
     - *reaction_prefix* [default='R _'] the prefix to add to input reaction ID's
     - *has_header* [default=False] if there is a header row in the csv file

    """

    assert os.path.exists(model_file), '\nPlease supply a file name that exists.'
    Rdict = getReactions(
        model_file,
        reaction_prefix=reaction_prefix,
        has_header=has_header,
        save_rpt=True,
    )
    R = parseReactions(Rdict)
    del Rdict
    if bounds_file != None and os.path.exists(bounds_file):
        B = getBounds(
            bounds_file, reaction_prefix=reaction_prefix, has_header=has_header
        )
        print('\nBounds loaded from file: {}\n'.format(bounds_file))
    else:
        B = {}

    addBoundsToReactions(R, B, default=BIG_BOUND)
    S = getSpecies(R)
    dumpReactionsToTxt(R, 'dump_reactions.txt')
    dumpSpeciesToTxt(S, 'dump_species.txt')

    fba = buildFBAobj(reactions=R, species=S, objfname=biomass_flux, modname=model_id)
    return fba


"""
if __name__ == '__main__':
    cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
    work_dir = os.path.join(cDir, 'spy')

    Bounds = getBounds(os.path.join(work_dir, 'Spy_bounds_semicolon.csv'))
    Reactions  = getReactions(os.path.join(work_dir, 'Spy_reaction_biomass_semicolon.csv'))
    ##  Bounds = getBounds(os.path.join(work_dir, 'boundaries_EFA.csv'))
    ##  Reactions  = getReactions(os.path.join(work_dir, 'reactions_EFA1.csv'))

    addBoundsToReactions(Reactions, Bounds, default=1000.0)
    dumpReactionsToTxt(Reactions, 'dump_reactions.txt')
    Species = getSpecies(Reactions)
    dumpSpeciesToTxt(Species, 'dump_species.txt')
    fba = buildFBAobj(Reactions, Species, 'biomass_LPL6_0', 'Spy')
"""
