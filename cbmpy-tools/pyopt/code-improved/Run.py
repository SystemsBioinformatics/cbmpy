"""
    pyOpt 1.0
    ==============
    
    A Python implementation of the OptKnock algorithm
    Copyright (C) 2016 Coco van Boxtel and Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands
    
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
    
    Authors: Coco van Boxtel and Brett G. Olivier
    Contact email: c.van.boxtel@vu.com
    Contact email: bgoli@users.sourceforge.net
    
"""

import os, sys
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pyOpt, productSelection, testIndependently

product_dict = {
'iaz_sbml3.xml':'R_EX_succ_e_',
'imm_sbml3.xml':'R_EX_succ_e_',
'toy_model_biofuel.l3.xml':'R09',
'iTM686_biofuels.sbml3.xml':'R_EX_etoh_e',
'ModelinRefCondition.xml':['R_EX_ac_e', 'R_EX_etoh_e']
}

biomass_dict = {
'iaz_sbml3.xml':'R_biomass_core',
'imm_sbml3.xml': 'R_biomass_SC5_notrace',
'toy_model_biofuel.l3.xml':'R17',
'iTM686_biofuels.sbml3.xml':'R_BiomassHetero',
'ModelinRefCondition.xml':'R_BiomassAuto'
}


### User input ######

#modelFile = 'ModelinRefCondition.xml'
# modelFile = 'iTM686_biofuels.sbml3.xml'
# modelFile = 'toy_model_biofuel.l3.xml'
# modelFile = 'iaz_sbml3.xml'
# modelFile = 'imm_sbml3.xml'


# percentage of bio_reaction that should be maintained
#objMinFactor = 0.1

# the minimum number of fluxes that need to be potentially active (total number of fluxes - deletions)
#maxDelete = 3
    
# Define value that indicates 'unbounded' in sbml file, to be translated to cplex.inf
# If boundaries are already defined infinite, make sure this value is high
infinityValue = 99999
    
# Use gene knockouts instead of reaction knockouts
# Warning, the model needs to have well defined GPR associations; '(G1 and G2 or G3) or (G6)'
# USE_GENE = False
    
# Reduce number of knockouts
# Keep alpha value small enough to avoid potential dilution of the bilevel objective
# Can only be used if bilevelObjective is maximization
USE_KNOCKOUT_WEIGHTING = False
KNOCKOUT_WEIGHTING_ALPHA = 0.0002
    
# Allow a near optimal solution (e.g. 5% less than maximum)
# This will increase the solution speed.
SOLUTION_FROM_OPTIMUM = 0.0

#####################

#bio_reaction = biomass_dict[modelFile]
#bilevelObjectives = product_dict[modelFile]


# prompt, uncomment this block if using User input
print 'Give name of file, located in the folder "Models"'
modelFile = raw_input('> ')

print 'What is the reaction id of the biomass forming reaction?'
bio_reaction = raw_input('> ')

print 'Do you automatically want to search for product candidates?'
autoProd = raw_input('> ')
if autoProd == 'Yes' or autoProd =='yes':
    bilevelObjectives, byproducts = productSelection.findProduct(modelFile, bio_reaction)
else:
    print 'Which reaction do you want to optimize for then?'
    bilevelObjectives = raw_input('> ')
    bilevelObjectives = [bilevelObjectives]


print '\nWhat fraction of biomass production do you want to maintain?'
objMinFactor = float(raw_input('> '))

print 'What is the maximum amount of knockouts you would like to have?'
maxDelete = float(raw_input('> '))

print 'Do you want to know gene knockouts?'
doGenes = raw_input('> ')
if doGenes == 'Yes' or doGenes =='yes':
    USE_GENE = True
else:
    USE_GENE = False



# Open new file for results
resultFile = file(os.path.join('Results', modelFile.replace('.xml', '') + '.result.txt'),'w')
resultFile.write('=================\n')
resultFile.write('pyOpt Solutions\n=================\n')
resultFile.write('\nMaximum deletions: {}'.format(maxDelete))
if len(bilevelObjectives) > 1:
    resultFile.write('\nAll metabolites tested: {}'.format(byproducts))
if USE_KNOCKOUT_WEIGHTING:
    resultFile.write('\nKnockout weighting (minimization): ON')
    resultFile.write('\nWARNING! The outer optimum value is slightly lower, due to use of knockout weighting')
resultFile.write('\n\n=================\n\n')


# Run pyOpt for compounds of interest
for reaction in bilevelObjectives:
    if reaction > 1:
        [bilevelObjective, val1, bio_reaction, val2, delGen,delReact] = pyOpt.runOptKnock('SynthSinks/'+reaction[1], (reaction[0],1), bio_reaction, objMinFactor, maxDelete, infinityValue, USE_GENE, USE_KNOCKOUT_WEIGHTING, KNOCKOUT_WEIGHTING_ALPHA, SOLUTION_FROM_OPTIMUM)
    else:
        [bilevelObjective, val1, bio_reaction, val2, delGen,delReact] = pyOpt.runOptKnock(modelFile, (reaction,1), bio_reaction, objMinFactor, maxDelete, infinityValue, USE_GENE, USE_KNOCKOUT_WEIGHTING, KNOCKOUT_WEIGHTING_ALPHA, SOLUTION_FROM_OPTIMUM)

    if abs(round(val1,5)) != 0.0:
        # test outcome independent from OptKnock implementation
        value1,value2,value3 = testIndependently.testKnockouts('SynthSinks/'+reaction[1], reaction[0], bio_reaction, USE_GENE,delGen,delReact)

        # write results in file
        resultFile.write('\nProduct reaction ({}): min {}, max {}'.format(bilevelObjective[0],round(value1,5),round(value2,5)))
        resultFile.write('\nBiomass reaction ({}): {}'.format(bio_reaction,round(value3,5)))
        if USE_GENE:
            resultFile.write('\nDeleted genes ({}): {}'.format(len(delGen), delGen))
        else:
            resultFile.write('\nDeleted reactions ({}): {}'.format(len(delReact), delReact))
        resultFile.write('\n_________________\n')

resultFile.close()







