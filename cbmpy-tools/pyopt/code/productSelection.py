import pyscescbm as cbm
import sys

byproducts = []

def getByProd(modelFile,coupledSubstrates, oldSubstrates, reactionList):
    newSubstrates = []
    # Find product associated reactions
    for product in coupledSubstrates:
        if not ('M_atp_' in product or 'M_adp_' in product or 'M_h2o' in product or 'M_nadph' in product or 'M_nadh' in product or 'M_nad' in product or 'M_h_' in product or 'M_pi' in product or 'M_ppi_' in product or 'M_no4_'in product or 'M_so3' in product or 'M_nh4' in product or 'M_na1' in product or 'M_Fe2' in product or 'M_fdxr' in product or 'M_fdxo' in product or 'M_ACP' in product or 'M_amp' in product or 'M_co2' in product or 'M_coa' in product or 'M_accoa' in product or 'M_glu' in product or 'M_pyr_' in product):
            reactions = modelFile.getFluxesAssociatedWithSpecies(product)
            for reaction in reactions:
                # Verify that this is an unique assay
                if reaction not in reactionList:
                    allProducts = modelFile.getReaction(reaction[0]).getProductIds()
                    # Select producing reaction
                    if product in allProducts:
                        reactionList.append(reaction[0])
                        # Define precursors and byproducts
                        subs = modelFile.getReaction(reaction[0]).getSubstrateIds()
                        for substrate in subs:
                                newSubstrates.append(substrate)
                        for extra in allProducts:
                            if not extra in oldSubstrates and extra != product:
                                if not ('M_atp_' in extra or 'M_adp_' in extra or 'M_h2o' in extra or 'M_nadph' in extra or 'M_nadh' in extra or 'M_nad' in extra or 'M_h_' in extra or 'M_pi' in extra or 'M_ppi_' in extra or 'M_no4_'in extra or 'M_so3' in extra or 'M_nh4' in extra or 'M_na1' in extra or 'M_Fe2' in extra or 'M_fdxr' in extra or 'M_fdxo' in extra or 'M_ACP' in extra or 'M_amp' in extra or 'M_co2' in extra or 'M_coa' in extra or 'M_accoa' in extra or 'M_glu' in extra or 'M_pyr_' in extra or 'M_co_' in extra) and 'C' in modelFile.getSpecies(extra).getChemFormula():
                                    byproducts.append(extra)
                                    print  product + ' with as byproduct\t' + str(modelFile.getSpeciesIds(extra))

                
    newSubstrates = list(set(newSubstrates)-set(oldSubstrates))
    oldSubstrates = oldSubstrates + newSubstrates

    return reactionList, byproducts, newSubstrates, oldSubstrates


def createSink(modelFile, to_sink):
    
    sinkR = modelFile.getReactionIds('R_EX')
    sinkM = {}
    sinkIDs = []
    
    # Define metabolites in cytoplasm with exchange reaction
    for exR in sinkR:
        reaction = modelFile.getReaction(exR)
        sinkM[reaction.getSpeciesIds()[0][:-2]+'_c'] = exR
    
    # Verify if exchange reaction excist for products of interest
    # Otherwise create sink reaction
    for byproduct in to_sink:
        if byproduct in sinkM:
            sinkIDs.append(sinkM[byproduct])
        else:
            try:
                s = modelFile.getSpecies(byproduct)
                speciesname = s.getName()
                sinkid = 'R_EX_' + byproduct[2:]
                sinkname = speciesname + '_sink'
                modelFile.createReaction(sinkid, name = sinkname, reversible=True, create_default_bounds=True)
                modelFile.createReactionReagent(sinkid, byproduct, -1)
                modelFile.setReactionBounds(sinkid,0,999999)
                modelFile.getReaction(sinkid).setAnnotation('SUBSYSTEM','Exchange reactions')
                modelFile.getReaction(sinkid).setAnnotation('NAME',sinkname)
                sinkIDs.append(sinkid)
            except:
                print 'No sink reaction created for: '+byproduct
    

    return modelFile, sinkIDs

def findProduct(modelFile, biomass_reaction):
    
    reactionList = []
    oldSubstrates = []
    
    fileName = modelFile[:-4]
    modelFile = cbm.CBRead.readSBML3FBC('Models/'+modelFile)
    newFileName = fileName + '_extraSinks.xml'

    # Take biomass compounds as first products to be analysed
    bioSubstrates = modelFile.getReaction(biomass_reaction).getSubstrateIds()
    reactionList, byproducts, newSubstrates, oldSubstrates = getByProd(modelFile, bioSubstrates, oldSubstrates, reactionList)

    # Repeat procedure for substrate metabolites found
    while len(newSubstrates) > 0:
        reactionList, byproducts, newSubstrates, oldSubstrates = getByProd(modelFile, newSubstrates, oldSubstrates, reactionList)

    modelFile, sinkIDs = createSink(modelFile,list(set(byproducts)))

    cbm.CBWrite.writeSBML3FBC(modelFile, newFileName, directory='Models')

    return newFileName, sinkIDs



