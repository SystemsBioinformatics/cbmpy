"""
Script for setting a growth condition to a sbml model

"""

import pyscescbm as cbm
import sys

def setCondition(modelRaw, bio_reaction):

    autotrophicGrowth = False
    mixotrophicGrowth = False
    heterotrophicGrowth = False

    if bio_reaction == 'R_BiomassAuto':
        autotrophicGrowth = True
    elif bio_reaction == 'R_BiomassMixo':
        mixotrophicGrowth = True
    elif bio_reaction == 'R_BiomassHetero':
        heterotrophicGrowth = True
    else:
        print '\nSetting growth condition failed. Check if the condition input variable is properly defined.\n'
        sys.exit()

    try:
        # O2 uptake and excretion are not limited

        # Carbon Limited state
        # Adjusted to settings model
        # Objective value without knockouts .0523306361225665
        if autotrophicGrowth:
            print 'Conditions adjusted to Autotrophic'
            modelRaw.setObjectiveFlux(bio_reaction, osense = 'maximize')
            modelRaw.setReactionBounds('R_EX_photon_LPAREN_e_RPAREN_', -30, 0)
            modelRaw.setReactionBounds('R_EX_hco3_LPAREN_e_RPAREN_', -3.7, 0)
            modelRaw.setReactionBounds('R_EX_glc_LPAREN_e_RPAREN_', 0, 99999)
           

        # Carbon Limited state
        if mixotrophicGrowth:
            print 'Conditions adjusted to Mixotrophic'
            modelRaw.setObjectiveFlux(bio_reaction, osense = 'maximize')
            modelRaw.setReactionBounds('R_EX_photon_LPAREN_e_RPAREN_', -30, 0)
            modelRaw.setReactionBounds('R_EX_hco3_LPAREN_e_RPAREN_', -3.7, 0)
            modelRaw.setReactionBounds('R_EX_glc_LPAREN_e_RPAREN_', -1, 99999)


        if heterotrophicGrowth:
            print 'Conditions adjusted to Heterotrophic'
            modelRaw.setObjectiveFlux(bio_reaction, osense = 'maximize')
            modelRaw.setReactionBounds('R_EX_photon_LPAREN_e_RPAREN_', 0, 0)
            modelRaw.setReactionBounds('R_EX_hco3_LPAREN_e_RPAREN_', 0, 0)
            modelRaw.setReactionBounds('R_EX_glc_LPAREN_e_RPAREN_', -1, 99999)

    except:
        print '\nSetting growth condition failed. Check the specifications of the growth condition in GrowthCondition.py\n'
        sys.exit()


    cbm.CBSolver.analyzeModel(modelRaw)

    return modelRaw

