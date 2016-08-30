import cbmpy as cbm
import pandas as pd
import sys

def printReactions(modID, onlyExch=False, suff=None):

    """
    Prints reaction bounds, reactions' stoichiometry and their values
    """

    for ri in modID.reactions:

        if onlyExch:
            if ri.getPid().startswith('R_EX_'):
                print modID.getReactionBounds(ri.getPid()), modID.getReaction(ri.getPid()).getStoichiometry(), \
                    modID.getReactionValues()[ri.getPid()]

        elif suff is not None:
            if ri.getPid().endswith(suff):
                print modID.getReactionBounds(ri.getPid()), modID.getReaction(ri.getPid()).getStoichiometry(), \
                    modID.getReactionValues()[ri.getPid()]

        else:
            print modID.getReactionBounds(ri.getPid()), modID.getReaction(ri.getPid()).getStoichiometry(), \
                modID.getReactionValues()[ri.getPid()]


def checkEndings(aList, suf1, suf2):

    # check whether there are strings in the list with the ending suf1 and suf2, respectively
    # at least one string of each type is required to return TRUE

    boolSuf1 = any([x.endswith(suf1) for x in aList])
    boolSuf2 = any([x.endswith(suf2) for x in aList])

    return all([boolSuf1, boolSuf2])

def getReaBySuffixesMod(fba_mod, suf1, suf2):

    reaSpecDict = {reai.getPid(): fba_mod.getReaction(reai.getPid()).getSpeciesIds() for reai in fba_mod.reactions}

    newDict = {k: checkEndings(v, suf1, suf2) for k, v in reaSpecDict.iteritems()}

    return filter(lambda x: newDict[x], newDict.iterkeys())

def get_lb(val, tol_bound=0.3):

    if val < 0:

        return val * (1. + tol_bound)

    return val * (1. - tol_bound)

def get_ub(val, tol_bound=0.3):

    if val < 0:

        return val * (1. - tol_bound)

    return val * (1. + tol_bound)


bulgTherOrg = cbm.readSBML3FBC('bulgTherm.xml')

# fixing model
bulgTherOrgFix = bulgTherOrg.clone()

# add glutamine export
bulgTherOrgFix.createReaction('R_gln_export_bulg', 'export of glutamine')
# define the reagents in the created reaction
bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_gln_L_c_bulg', 1)
bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_h_c_bulg', 1)
bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_gln_L_e', -1)
bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_h_e', -1)
bulgTherOrgFix.setReactionBounds('R_gln_export_bulg', cbm.NINF, 0)
## # bgoli:20160307 turned reaction around so that it is a positive export (re-signed coiefficients and swapped bounds)
## bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_gln_L_c_bulg', -1)
## bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_h_c_bulg', -1)
## bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_gln_L_e', 1)
## bulgTherOrgFix.createReactionReagent('R_gln_export_bulg', 'M_h_e', 1)
## bulgTherOrgFix.setReactionBounds('R_gln_export_bulg', 0, cbm.INF)

# get all metabolites with an L in their name (should be mainly amino acids)
allL = [idi for idi in bulgTherOrgFix.getSpeciesIds() if '_L_c_' in idi]

# these are not amino acids and we don't want to add source reactions for them
ign = ['M_2hxic_L_c_bulg',
       'M_2hxic_L_c_ther',
       'M_citr_L_c_ther',
       'M_cysth_L_c_bulg',
       'M_cysth_L_c_ther',
       'M_hcys_L_c_bulg',
       'M_hcys_L_c_ther',
       'M_hom_L_c_bulg',
       'M_hom_L_c_ther',
       'M_lac_L_c_bulg',
       'M_lac_L_c_ther',
       'M_mal_L_c_bulg',
       'M_mal_L_c_ther',
       'M_orn_L_c_ther',
       'M_pser_L_c_ther',
       'M_scys_L_c_ther'
       ]

# get all amino acids
aaAdd = [aai for aai in allL if aai not in ign]
aaAdd.append('M_gly_c_bulg')
aaAdd.append('M_gly_c_ther')

# add a source reaction for each amino acid
for aai in aaAdd:
    if '_bulg' in aai:
        cbm.CBTools.addSourceReaction(bulgTherOrgFix, aai, lb=0., ub=cbm.INF)

bulgTherOrgFix.changeAllFluxBoundsWithValue(999999.0, cbm.INF)
bulgTherOrgFix.changeAllFluxBoundsWithValue(-999999.0, cbm.NINF)

bulgTherOrgFix.setReactionBounds('R_EX_cit_e_', 0, cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_fol_e_', 0., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_ac_e_', 0., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_ura_e_', 0., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_urea_e_', 0., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_pydam_e_', -1000, cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_thm_e_', -1000, cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_ribflv_e_', -1000, cbm.INF)

bulgTherOrgFix.setReactionBounds('R_EX_xan_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_gua_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_hxan_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_pnto_R_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_nac_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_btn_e_', -1000., cbm.INF)
bulgTherOrgFix.setReactionBounds('R_EX_cas', 0, 0)

cbm.writeSBML3FBC(bulgTherOrgFix, 'bulgTherSrcOnlyBulg.xml', directory='models')

sys.exit(0)

# get all reactions between cytosol and external compartment
bulgReaEx = getReaBySuffixesMod(bulgTherOrgFix, '_e', '_c' + '_bulg')
therReaEx = getReaBySuffixesMod(bulgTherOrgFix, '_e', '_c' + '_ther')

# filter all reactions between _c and _e in which an amino acid participates
aaReaMinSet = [reai for reai in bulgReaEx + therReaEx if
         any([si in aaAdd for si in bulgTherOrgFix.getReaction(reai).getSpeciesIds()])]

minSet = pd.DataFrame({'minRea': aaReaMinSet})

# sys.exit(0)

minSet.to_csv('minimizationsets/minset.csv', index=False, header=False)

# get the constrains
constrDF = pd.read_csv('constraints_EV.csv')

metReaAss = {}
metReaAss['Lactate_(mM)'] = {'reaID': 'R_EX_lac_L_e_', 'metID': "M_lac_L_e"}
metReaAss['Valine_(uM)'] = {'reaID': 'R_EX_val_L_e_', 'metID': "M_val_L_e"}
metReaAss['Galactose_(mM)'] = {'reaID': "R_EX_gal_e_", 'metID': "M_gal_e"}
metReaAss['Serine_(uM)'] = {'reaID': 'R_EX_ser_L_e_', 'metID': 'M_ser_L_e'}
metReaAss['Glutamine_(uM)'] = {'reaID': 'R_EX_gln_L_e_', 'metID': 'M_gln_L_e'}
metReaAss['Methionine_(uM)'] = {'reaID': 'R_EX_met_L_e_', 'metID': 'M_met_L_e'}
metReaAss['Tryptophan_(uM)'] = {'reaID': 'R_EX_trp_L_e_', 'metID': 'M_trp_L_e'}
metReaAss['Phenylalanine_(uM)'] = {'reaID': 'R_EX_phe_L_e_', 'metID': 'M_phe_L_e'}
metReaAss['Isoleucine_(uM)'] = {'reaID': 'R_EX_ile_L_e_', 'metID': 'M_ile_L_e'}
metReaAss['Lysine_(uM)'] = {'reaID': 'R_EX_lys_L_e_', 'metID': 'M_lys_L_e'}
metReaAss['Leucine_(uM)'] = {'reaID': 'R_EX_leu_L_e_', 'metID': 'M_leu_L_e'}
metReaAss['Lactose_(mM)'] = {'reaID': 'R_EX_lcts_e_', 'metID': 'M_lcts_e'}
metReaAss['sumQpcrBul'] = {'reaID': "R_biomass_LBUL10_bulg", 'metID': None}
metReaAss['sumQpcrTher'] = {'reaID': "R_biomass_STR_ther", 'metID': None}
metReaAss['Pyruvate_(mM)'] = {'reaID': "R_EX_pyr_e_", 'metID': 'M_pyr_e'}  # exchange reaction does not exist
metReaAss['Succinate_(mM)'] = {'reaID': "R_EX_succ_e_", 'metID': 'M_succ_e'}
metReaAss['Formate_(mM)'] = {'reaID': "R_EX_for_e_", 'metID': 'M_for_e'}
metReaAss['Histidine_(uM)'] = {'reaID': "R_EX_hist_L_e_", 'metID': 'M_hist_L_e'}
metReaAss['Arginine_(uM)'] = {'reaID': "R_EX_arg_L__e_", 'metID': 'M_arg_L_e'}
metReaAss['Alanine_(uM)'] = {'reaID': "R_EX_ala_L_e_", 'metID': 'M_ala_L_e'}
metReaAss['Glycine_(uM)'] = {'reaID': "R_EX_gly_e_", 'metID': 'M_gly_e'}
metReaAss['Acetate_(mM)'] = {'reaID': "R_EX_ac_e_", 'metID': 'M_ac_e'}
metReaAss['sumQpcrBulLog'] = {'reaID': "R_biomass_LBUL10_bulg_log", 'metID': None}
metReaAss['sumQpcrTherLog'] = {'reaID': "R_biomass_STR_ther_log", 'metID': None}

dfCons = constrDF.iloc[60:62:200, :]
# drop the columns that do not represent fluxes and won't be constraint
dfCons.drop(['Time_[h]', 'sumQpcrBul', 'sumQpcrTher', 'R_biomass_LBUL10_bulg_log', 'R_biomass_STR_ther_log'],
            1, inplace=True)

# drop the weirdly fitted metabolites
dfCons.drop([metReaAss['Pyruvate_(mM)']['reaID'], metReaAss['Acetate_(mM)']['reaID']], 1, inplace=True)

dfCons = dfCons.T

dfCons.rename(columns={dfCons.columns[0]: 'actualCons'}, inplace=True)


dfCons = dfCons.merge(dfCons.actualCons.apply(lambda s: pd.Series({'lb': get_lb(s), 'ub': get_ub(s)})),
    left_index=True, right_index=True)

dfCons.drop('actualCons', 1, inplace=True)

# avoid numerical issues caused by too narrow bounds
dfCons = dfCons*1000

dfCons.to_csv('constraints/consets.csv', index=True, header=False)

sys.exit(0)

rea_and_bounds = [list(bulgTherOrgFix.getReactionBounds(reai)[0:-1]) for reai in bulgReaEx + therReaEx]

dfReaBounds = pd.DataFrame(rea_and_bounds, columns=['reaID', 'lb', 'ub'])
dfReaBounds.set_index('reaID', inplace=True)
dfReaBounds = pd.concat([dfCons, dfReaBounds])
dfReaBounds.to_csv('constraints/consets.csv', index=True, header=False)

