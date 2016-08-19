import pyscescbm as cbm

# Set values
biomass = 'R_BiomassAuto'
fuel = 'R_EX_ade_c'

knockouts = ['sll0554', 'sll1430']
# Set growth conditions
model = cbm.CBRead.readSBML3FBC('Models/ModelinRefCondition_extraSinks.xml')
model.splitEqualityFluxBounds()
model.createGeneAssociationsFromAnnotations()


# Apply Knockouts
for g_ in knockouts:
        #model.setReactionBounds(g_, 0.0, 0.0)
        model.setGeneInactive(g_, update_reactions=True)

cbm.CBSolver.analyzeModel(model)

ObjValue = model.getReaction(biomass).getValue()

model.setReactionLowerBound(biomass, ObjValue)
model.setReactionUpperBound(biomass, ObjValue)

cbm.CBSolver.analyzeModel(model, oldlpgen=False)

test = cbm.FluxVariabilityAnalysis(model, selected_reactions = [fuel])
print 'Biomass: ' + str(ObjValue)
print 'Minimum: ' + str(test[0][0][2])
print 'Maximum: ' + str(test[0][0][3])

