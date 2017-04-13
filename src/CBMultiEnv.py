"""
CBMPy: CBMultiEnv module
========================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2017 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

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
Last edit: $Author: bgoli $ ($Id: CBMultiEnv.py 575 2017-04-13 12:18:44Z bgoli $)

"""

# preparing for Python 3 port
from __future__ import division, print_function
from __future__ import absolute_import
#from __future__ import unicode_literals

import os, time, numpy
from . import CBRead, CBWrite, CBTools, CBNetDB, CBSolver, CBMultiCore

HAVE_XLWT = False
try:
    import xlwt
    HAVE_XLWT = True
except:
    pass

from .CBConfig import __CBCONFIG__ as __CBCONFIG__
__DEBUG__ = __CBCONFIG__['DEBUG']
__version__ = __CBCONFIG__['VERSION']

#print('\n***\nMultiple Environment Module ({})\n***\n'.format(__version__))

def addProteinCostAnnotation(F, reaction_gene_map, gene_peptide_map, gene_annotation_key='GENE ASSOCIATION'):
    reaction_gene_map = reaction_gene_map.copy()
    gene_peptide_map = gene_peptide_map.copy()
    for r in F.reactions:
        if reaction_gene_map[r.getId()] != None:
            r.annotation.update({'CBM_PEPTIDE_LENGTHS' : [(g, gene_peptide_map[g]) for g in reaction_gene_map[r.getId()]]})
        else:
            r.annotation.update({'CBM_PEPTIDE_LENGTHS' : []})
        GA = None
        if gene_annotation_key in r.annotation:
            GA = gene_annotation_key
        elif 'GENE ASSOCIATION' in r.annotation:
            GA = 'GENE ASSOCIATION'
        elif 'GENE_ASSOCIATION' in r.annotation:
            GA = 'GENE_ASSOCIATION'
        if GA != None:
            newAss = r.annotation[GA]
            for rl in r.annotation['CBM_PEPTIDE_LENGTHS']:
                ##  print rl
                ##  print newAss
                newAss = newAss.replace(rl[0], str(rl[1]))
            newAss = newAss.replace('and','+').replace('(','').replace(')','')
            r.annotation.update({'CBM_PEPTIDE_SUBSTITUTED_LENGTHS' : newAss})
            ##  print r.annotation
            if 'or' in r.annotation[GA]:
                ##  print ' '
                ##  print r.annotation[GA]
                ##  print r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS']
                ##  print r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'].split('or')
                splitor = []
                for ass in r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'].split('or'):
                    exec('splitor.append(%s)' % ass)
                ##  print splitor
                r.annotation.update({'CBM_PEPTIDE_LENGTH_MIN' : min(splitor)})
                r.annotation.update({'CBM_PEPTIDE_LENGTH_MAX' : max(splitor)})
            else:
                if r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'] != None and r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'] != '':
                    maxL = 0
                    ##  print r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS']
                    print(r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'])
                    try:
                        exec('maxL = %s' % r.annotation['CBM_PEPTIDE_SUBSTITUTED_LENGTHS'])
                    except:
                        time.sleep(1)
                    r.annotation.update({'CBM_PEPTIDE_LENGTH_MIN' : maxL})
                    r.annotation.update({'CBM_PEPTIDE_LENGTH_MAX' : maxL})
                else:
                    r.annotation.update({'CBM_PEPTIDE_LENGTH_MIN' : None})
                    r.annotation.update({'CBM_PEPTIDE_LENGTH_MAX' : None})
            ##  print r.getId(), r.annotation['CBM_PEPTIDE_LENGTH_MIN'], r.annotation['CBM_PEPTIDE_LENGTH_MAX']
        else:
            r.annotation.update({'CBM_PEPTIDE_SUBSTITUTED_LENGTHS' : None})
            r.annotation.update({'CBM_PEPTIDE_LENGTH_MIN' : None})
            r.annotation.update({'CBM_PEPTIDE_LENGTH_MAX' : None})

def scaleProteinCostAnnotation(fba):
    LEN_MAX = max([r.annotation['CBM_PEPTIDE_LENGTH_MAX'] for r in fba.reactions if r.annotation['CBM_PEPTIDE_LENGTH_MAX'] != None])
    LEN_AVG = numpy.average([r.annotation['CBM_PEPTIDE_LENGTH_MAX'] for r in fba.reactions if r.annotation['CBM_PEPTIDE_LENGTH_MAX'] != None])

    # cost function (AAlen/avg(AAlen) or 1.0
    for R in fba.reactions:
        if R.annotation['CBM_PEPTIDE_LENGTH_MAX'] == None:
            R.annotation['CBM_PEPTIDE_COST'] = 1.0
        else:
            R.annotation['CBM_PEPTIDE_COST'] = R.annotation['CBM_PEPTIDE_LENGTH_MAX']/LEN_AVG
        R.annotation['CBM_AVG_PEPTIDE_LENGTH'] = LEN_AVG
    print('\nAverage protein length:', LEN_AVG)
        ##  print R.getId(), R.annotation['CBM_PEPTIDE_COST'], R.annotation['CBM_PEPTIDE_LENGTH_MAX']

# Standard amino acid abbreviations with ATP formation/production cost: Stouthamer (1973) Table 4 page 555
amino_acid_table = [['A', 'Ala', 'Alanine', '+1'],
                    ['R', 'Arg', 'Arginine', '-3'],
                    ['N', 'Asn', 'Asparagine', '-2'],
                    ['D', 'Asp', 'Aspartic acid', '+0'],
                    ['C', 'Cys', 'Cysteine', '-3'],
                    ['E', 'Glu', 'Glutamic acid', '+1'],
                    ['Q', 'Gln', 'Glutamine', '+0'],
                    ['G', 'Gly', 'Glycine', '+0'],
                    ['H', 'His', 'Histidine', '-7'],
                    ['I', 'Ile', 'Isoleucine', '-1'],
                    ['L', 'Leu', 'Leucine', '+3'],
                    ['K', 'Lys', 'Lysine', '+0'],
                    ['M', 'Met', 'Methionine', '-4'],
                    ['F', 'Phe', 'Phenylalanine', '-2'],
                    ['P', 'Pro', 'Proline', '+0'],
                    ['S', 'Ser', 'Serine', '+0'],
                    ['T', 'Thr', 'Threonine', '-2'],
                    ['W', 'Trp', 'Tryptophan', '-5'],
                    ['Y', 'Tyr', 'Tyrosine', '-2'],
                    ['V', 'Val', 'Valine', '+2']
                    ]


def runMultiStateFBA(input_list, dev_factor, options, dofva=False, fva_selected_reactions=None):
    if options['algorithm'] == "ABS_L1_METHOD":
        if options['with_cost']:
            work_dir = os.path.join(options['base_work_dir'], 'AbsL1_cost')
        else:
            work_dir = os.path.join(options['base_work_dir'], 'AbsL1_nocost')
    elif options['algorithm'] == "L1_METHOD":
        if options['with_cost']:
            work_dir = os.path.join(options['base_work_dir'], 'L1_cost')
        else:
            work_dir = os.path.join(options['base_work_dir'], 'L1_nocost')

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    IO_list = []
    ModelPrefixStrIO = ''
    for I in input_list:
        inputB = I[0]
        prefix = I[1]
        environ = None
        try:
            environ = I[2]
            print('INFO: custom environment defined:', environ)
        except:
            print('INFO: no custom environment defined')
        ModelPrefixStrIO += prefix
        ModelName = '%s(%s)' % (options['model_name_base'], inputB)

        fmg1 = CBRead.readSBML2FBA(options['input_model_file'], options['model_dir'])
        # get active (1st) objfunc
        biomass_name = fmg1.objectives[0].fluxObjectives[0].reaction
        # set name and inputs
        fmg1.id = prefix+options['model_name_base']
        # set environment
        if environ != None:
            environ = tuple(environ)
            for j_ in environ:
                print('Reaction', j_[0])
                print('LB', j_[1])
                print('UB', j_[2])
                assert fmg1.getReaction(j_[0]) != None, '\n%s is not a valid reaction id!' % j_[0]
                fmg1.setReactionLowerBound(j_[0], j_[1])
                fmg1.setReactionUpperBound(j_[0], j_[2])
            #os.sys.exit()
        # prefix names
        fmg1.setPrefix(prefix, 'all')
        inputB = prefix+inputB
        # write LP
        CBTools.addStoichToFBAModel(fmg1)

        #cpx1 = CBSolver.cplx_constructLPfromFBA(fmg1, fname=None)
        LPFile1 = CBWrite.writeModelLP(fmg1, work_dir)
        cpx1 = CBSolver.cplx_getModelFromLP(LPFile1)
        CBSolver.cplx_setObjective(cpx1, prefix+'obj1', [(1.0, inputB)], 'max', reset=True)
        CBSolver.cplx_setBounds(cpx1, inputB, min=-options['myInf'], max=0.0)
        CBSolver.cplx_setSingleConstraint(cpx1, 'C_'+prefix+biomass_name, expr=[(1, prefix+biomass_name)], sense='E', rhs=options['target_biomass'])
        CBSolver.cplx_writeLPtoLPTfile(cpx1, os.path.join(work_dir, prefix+'%s(%s)_Imin' % (options['model_name_base'], I[0])), title='%s input minimization' % I[0], Dir=None)
        CBSolver.cplx_Solve(cpx1)
        fba_sol1, objf_name1, objf_val1 = CBSolver.cplx_getOptimalSolution(cpx1)
        ##  for j in fba_sol:
            ##  print '%s: %s' % (j, fba_sol[j])
        print('\nCPLEX returns objective {} = {}'.format(objf_name1, objf_val1))
        slv_stat = CBSolver.cplx_getSolutionStatus(cpx1)
        print('\nStatus: {}'.format(slv_stat))
        IO_list.append((I[0], objf_val1, I[1], slv_stat, environ))

    print('\nInput values for target: {} = {}'.format(biomass_name, options['target_biomass']))
    F = file(os.path.join(work_dir, 'find_input(%s)_(%s).csv' % (ModelPrefixStrIO, options['target_biomass'])), 'w')
    F.write('input,value,prefix,status,environ\n')
    for IO in IO_list:
        CBTools.processBiGGchemFormula(fmg1)
        CBTools.processSBMLAnnotationNotes(fmg1, annotation_key='note')

        ##  print '%s, %.20f, %s, %s\n' % IO
        F.write('%s,%3.8f,%s,%s,\"%s\"\n' % IO)

    F.flush()
    F.close()

    DB = None
    if options['with_cost']:
        ##  db_cursor, db_conn = CBTools.ConnectDB('ecoli.sqlite', work_dir=options['DB_dir'])
        DB = CBNetDB.DBTools()
        DB.connectSQLiteDB('ecoli.sqlite', options['DB_dir'])

    multisolve = []
    for IO in IO_list:
        print(IO[3])
        assert IO[3] == 'LPS_OPT', '\nBAD LP ...'
        inputB = IO[0]
        inputValue = IO[1]
        environ = None
        try:
            print('Environ', IO[4])
            environ = IO[4]
            print('INFO: custom environment defined:', environ)
        except:
            print('INFO: no custom environment defined')
        ##  inputValue = round(IO[1], 8) # note this is a thumbsuck
        prefix = IO[2]
        fmg2 = CBRead.readSBML2FBA(options['input_model_file'], work_dir=options['model_dir'])
        # fix annotations
        CBTools.processBiGGchemFormula(fmg2)
        CBTools.processSBMLAnnotationNotes(fmg2, annotation_key='note')


        if options['with_cost']:
            # process GENE ANNOTATIONS USE initial version of model
            react_gene = CBTools.getModelGenesPerReaction(fmg2, gene_pattern=options['gene_pattern'])
            gene_react, ecoli_genes, no_gene = CBTools.getReactionsPerGene(react_gene)
            # recreate KeGG identifiers
            ecoli_genes = [options['model_KeGG_prefix']+g for g in ecoli_genes]
            # retrieve and process amino acid sequences

            # extract info from database
            gene_peplen = {}
            for G in ecoli_genes:
                Glen = DB.fetchAll('SELECT aa_len FROM ecoli_genes WHERE gene="%s"' % G)
                if len(Glen) > 0 and len(Glen[0]) > 0:
                    Glen = Glen[0][0]
                    gene_peplen.update({G.replace(options['model_KeGG_prefix'],'') : Glen})
                else:
                    print('WARNING: no DB data for gene: {}'.format(G))
            ##  gene_peplen = CBTools.getPeptideLengthsFromDB(db_cursor, ecoli_genes, options['model_KeGG_prefix'])
            addProteinCostAnnotation(fmg2, react_gene, gene_peplen)
            scaleProteinCostAnnotation(fmg2)

        fmg2.setReactionBound(inputB, inputValue, 'lower')
        if options['relax_input_bound']:
            fmg2.setReactionBound(inputB, 0.0, 'upper')
        else:
            fmg2.setReactionBound(inputB, inputValue, 'upper')

        # set environment
        if environ != None:
            for j_ in environ:
                print('Reaction', j_[0])
                print('LB', j_[1])
                print('UB', j_[2])
                assert fmg2.getReaction(j_[0]) != None, '\n%s is not a valid reaction id!' % j_[0]
                fmg2.setReactionLowerBound(j_[0], j_[1])
                fmg2.setReactionUpperBound(j_[0], j_[2])
            #raw_input('Im bored. What play thing can you offer me today?')

        # prefix names
        fmg2.id = prefix+options['model_name_base']
        fmg2.setPrefix(prefix, 'all')

        CBTools.addStoichToFBAModel(fmg2)

        fmg2lp = CBSolver.cplx_constructLPfromFBA(fmg2, fname=None)
        #LPFile2 = CBWrite.writeModelLP(fmg2, work_dir=work_dir)
        #fmg2lp = CBSolver.cplx_getModelFromLP(LPFile2)
        CBSolver.cplx_Solve(fmg2lp)
        fmg2.objectives[0].solution, fmg2id, fmg2.objectives[0].value = CBSolver.cplx_getOptimalSolution(fmg2lp)
        slv_stat = CBSolver.cplx_getSolutionStatus(fmg2lp)
        print('\nStatus: {}'.format(slv_stat))
        print("{} = {}".format(fmg2id, fmg2.objectives[0].value))
        multisolve.append(fmg2)
        time.sleep(3)

    # quick and dirty check to see that the reactions are in the same order
    if options['run_order_test']:
        for M in range(1,len(multisolve)):
            for j in range(len(multisolve[0].N.col)):
                ##  print multisolve[0].N.col[j][2:], multisolve[M].N.col[j][2:]
                assert multisolve[0].N.col[j][2:] == multisolve[M].N.col[j][2:], '\nN col mismatch: %s not e %s' % (multisolve[0].N.col[j][2:], multisolve[M].N.col[j][2:])
            for j in range(len(multisolve[0].N.row)):
                ##  print multisolve[0].N.row[j][2:], multisolve[M].N.row[j][2:]
                assert multisolve[0].N.row[j][2:] == multisolve[M].N.row[j][2:], '\nN row mismatch: %s != %s' % (multisolve[0].N.row[j][2:], multisolve[M].N.row[j][2:])
            idx = multisolve[0].getReactionIds()
            test_idx = multisolve[M].getReactionIds()
            for j in range(len(idx)):
                print(idx[j][2:], test_idx[j][2:])
                assert idx[j][2:] == test_idx[j][2:], '\nReaction ID mismatch: %s != %s' % (idx[j][2:], test_idx[j][2:])
        print('\nOrder test completed succesfully')
        os.sys.exit()

    # get model prefixes in various forms for output happiness
    ModelPrefixes = [m.prefix for m in multisolve]
    ModelPrefixStr = ''
    for p in ModelPrefixes:
        ModelPrefixStr += p
    print(ModelPrefixStr)

    objective_base = tuple([f.objectives[0].value for f in multisolve])

    RESULT_DATA = []
    RESULT_LPS = []

    ##  RESULT = numpy.zeros((len(GENERIC_REACTION_OBJECTS), len(multisolve)*len(dev_factor)+1), numpy.double)
    for devF in range(len(dev_factor)):
        for f in range(len(multisolve)):
            multisolve[f].objectives[0].value = objective_base[f]*dev_factor[devF]
        if options['with_cost']:
            if options['algorithm'] == "ABS_L1_METHOD":
                ##  raise NotImplementedError
                MultiFNbase = '%s_P(%s)_T(%s)_ABSL1_C' % (options['input_model_file'], ModelPrefixStr, dev_factor[devF])
                LPFileG = CBWrite.writeMinDistanceLP_absL1(MultiFNbase, multisolve, work_dir=work_dir, with_protein_cost=True, ignoreDistance=[], bigM=options['L1_BIGM_VALUE'])
            elif options['algorithm'] == "L1_METHOD":
                MultiFNbase = '%s_P(%s)_T(%s)_L1_C' % (options['input_model_file'], ModelPrefixStr, dev_factor[devF])
                LPFileG = CBWrite.writeMinDistanceLP(MultiFNbase, multisolve, work_dir=work_dir, with_protein_cost=True)
        else:
            if options['algorithm'] == "ABS_L1_METHOD":
                MultiFNbase = '%s_P(%s)_T(%s)_ABSL1' % (options['input_model_file'], ModelPrefixStr, dev_factor[devF])
                LPFileG = CBWrite.writeMinDistanceLP_absL1(MultiFNbase, multisolve, work_dir=work_dir, with_protein_cost=False, ignoreDistance=[], bigM=options['L1_BIGM_VALUE'])
            elif options['algorithm'] == "L1_METHOD":
                MultiFNbase = '%s_P(%s)_T(%s)_L1' % (options['input_model_file'], ModelPrefixStr, dev_factor[devF])
                LPFileG = CBWrite.writeMinDistanceLP(MultiFNbase, multisolve, work_dir=work_dir, with_protein_cost=False)

        mGlp = CBSolver.cplx_getModelFromLP(LPFileG)
        TIME_START = time.time()
        if options['algorithm'] == "ABS_L1_METHOD":
            #mGlp.parameters.threads.set(4)
            CBSolver.cplx_SolveMILP(mGlp, auto_mipgap=True)
        elif options['algorithm'] == "L1_METHOD":
            CBSolver.cplx_Solve(mGlp)
        TIME_END = time.time()

        print('\nSanity check')
        for f in multisolve:
            print('Model: {}'.format(f.id))
            for inp in input_list:
                print(' {} lower: {}'.format(f.prefix+inp[0], f.getFluxBoundByReactionID(f.prefix+inp[0], 'lower').value))
                print(' {} upper: {}'.format(f.prefix+inp[0], f.getFluxBoundByReactionID(f.prefix+inp[0], 'upper').value))

        print('\nFBA optima')
        for i in multisolve:
            print(' {} has optimum: {}'.format(i.id, i.objectives[0].value))

        mGsol, mGobid, mGobjval = CBSolver.cplx_getOptimalSolution(mGlp)

        print('\n*****\nMultipleInput(%s states) Result'.format(len(multisolve)))
        print(" optimum: %s = %s".format(mGobid, mGobjval))
        print(' status: %s'.format(CBSolver.cplx_getSolutionStatus(mGlp)))
        print('*****')

        print('\nMultipleInput({} states) Solve Time'.format(len(multisolve)))
        print('Time: {}\nTime(m): {}\nTime(h): {}\n '.format(TIME_END-TIME_START, (TIME_END-TIME_START)/60.0, (TIME_END-TIME_START)/60.0**2))

        GENERIC_REACTION_OBJECTS = multisolve[0].reactions
        GENERIC_PREFIX = multisolve[0].prefix
        try:
            FlxAO_J = numpy.zeros((len(GENERIC_REACTION_OBJECTS), len(multisolve)), numpy.double)
            FlxAO_scaled = numpy.zeros((len(GENERIC_REACTION_OBJECTS), len(multisolve)+1), numpy.double)
        except AttributeError:
            FlxAO_J = numpy.zeros((len(GENERIC_REACTION_OBJECTS), len(multisolve)))
            FlxAO_scaled = numpy.zeros((len(GENERIC_REACTION_OBJECTS), len(multisolve)+1))
        FluxNamesOut = []
        FluxScaling = []

        for flx in range(len(GENERIC_REACTION_OBJECTS)):
            Rid = GENERIC_REACTION_OBJECTS[flx].id[len(GENERIC_PREFIX):]
            FluxNamesOut.append(Rid)
            if options['with_cost'] and 'CBM_PEPTIDE_COST' in GENERIC_REACTION_OBJECTS[flx].annotation:
                cost = GENERIC_REACTION_OBJECTS[flx].annotation['CBM_PEPTIDE_COST']
            else:
                cost = 1.0
            FluxScaling.append(cost)
            for p in range(len(ModelPrefixes)):
                FlxAO_J[flx, p] = mGsol[ModelPrefixes[p] + Rid]
                FlxAO_scaled[flx, p] = mGsol[ModelPrefixes[p] + Rid]*cost
            #print Rid, FlxAO_scaled[flx]
        ##  print FlxAO_scaled.shape, ModelPrefixes
        ##  print FlxAO_scaled[:3]


        # TODO: AUTOMATE !!!!!
        # FLUX DIFFERENCE SUM
        all_flux_diff_xxx = numpy.zeros(3)
        if len(multisolve) == 2:
            all_flux_diff = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,1])
            all_flux_diff_xxx = numpy.abs(numpy.abs(FlxAO_scaled[:,0]) - numpy.abs(FlxAO_scaled[:,1]))
        elif len(multisolve) == 3:
            all_flux_diff1 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,1])
            all_flux_diff2 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,2])
            all_flux_diff3 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,2])
            all_flux_diff = all_flux_diff1  + all_flux_diff2  + all_flux_diff3
            all_flux_diff_xxx1 = numpy.abs(numpy.abs(FlxAO_scaled[:,0]) - numpy.abs(FlxAO_scaled[:,1]))
            all_flux_diff_xxx2 = numpy.abs(numpy.abs(FlxAO_scaled[:,0]) - numpy.abs(FlxAO_scaled[:,2]))
            all_flux_diff_xxx3 = numpy.abs(numpy.abs(FlxAO_scaled[:,1]) - numpy.abs(FlxAO_scaled[:,2]))
            all_flux_diff_xxx = all_flux_diff_xxx1  + all_flux_diff_xxx2  + all_flux_diff_xxx3
        elif len(multisolve) == 4:
            all_flux_diff1 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,1])
            all_flux_diff2 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,2])
            all_flux_diff3 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,3])
            all_flux_diff4 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,2])
            all_flux_diff5 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,3])
            all_flux_diff6 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,3])
            all_flux_diff = all_flux_diff1  + all_flux_diff2  + all_flux_diff3 +\
                all_flux_diff4  + all_flux_diff5  + all_flux_diff6
        elif len(multisolve) == 5:
            all_flux_diff1 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,1])
            all_flux_diff2 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,2])
            all_flux_diff3 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,3])
            all_flux_diff4 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,4])
            all_flux_diff5 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,2])
            all_flux_diff6 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,3])
            all_flux_diff7 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,4])
            all_flux_diff8 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,3])
            all_flux_diff9 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,4])
            all_flux_diff10 = numpy.abs(FlxAO_scaled[:,3] - FlxAO_scaled[:,4])
            all_flux_diff = all_flux_diff1  + all_flux_diff2  + all_flux_diff3 + all_flux_diff4  + all_flux_diff5 +\
                all_flux_diff6 + all_flux_diff7  + all_flux_diff8  + all_flux_diff9 + all_flux_diff10
        elif len(multisolve) == 6:
            all_flux_diff1 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,1])
            all_flux_diff2 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,2])
            all_flux_diff3 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,3])
            all_flux_diff4 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,4])
            all_flux_diff5 = numpy.abs(FlxAO_scaled[:,0] - FlxAO_scaled[:,5])
            all_flux_diff6 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,2])
            all_flux_diff7 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,3])
            all_flux_diff8 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,4])
            all_flux_diff9 = numpy.abs(FlxAO_scaled[:,1] - FlxAO_scaled[:,5])
            all_flux_diff10 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,3])
            all_flux_diff11 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,4])
            all_flux_diff12 = numpy.abs(FlxAO_scaled[:,2] - FlxAO_scaled[:,5])
            all_flux_diff13 = numpy.abs(FlxAO_scaled[:,3] - FlxAO_scaled[:,4])
            all_flux_diff14 = numpy.abs(FlxAO_scaled[:,3] - FlxAO_scaled[:,5])
            all_flux_diff15 = numpy.abs(FlxAO_scaled[:,4] - FlxAO_scaled[:,5])
            all_flux_diff = all_flux_diff1  + all_flux_diff2  + all_flux_diff3 + all_flux_diff4  + all_flux_diff5 +\
                all_flux_diff6 + all_flux_diff7  + all_flux_diff8  + all_flux_diff9 + all_flux_diff10 + all_flux_diff11 + all_flux_diff12  + all_flux_diff13  + all_flux_diff14 + all_flux_diff15
        else:
            all_flux_diff = FlxAO_scaled[:,-1]

        FlxAO_scaled[:,-1] = all_flux_diff

        print('Flux difference', all_flux_diff.sum(), all_flux_diff_xxx.sum())

        if options['algorithm'] == "ABS_L1_METHOD":
            SSAall = ({},{},{})
        elif options['algorithm'] == "L1_METHOD":
            SSAall = CBSolver.cplx_getSensitivities(mGlp)
            CBWrite.writeSensitivitiesToCSV(SSAall, os.path.join(work_dir, '%s.ssa.csv' % MultiFNbase))

        DD = {'dev_factor' : dev_factor[devF],
              'model_file' : MultiFNbase,
              'solution_status' : CBSolver.cplx_getSolutionStatus(mGlp),
              'obj_val' : mGobjval,
              'obj_id' : mGobid,
              'solution_cplex' :  mGsol.copy(),
              'solution' : FlxAO_J,
              'solution_scaled' : FlxAO_scaled,
              'flux_names' : tuple(FluxNamesOut),
              'flux_diff' : all_flux_diff.sum(),
              'flux_diff_new': all_flux_diff_xxx.sum(),
              'flux_scaling' : tuple(FluxScaling),
              'sensitivity' : SSAall,
              'prefixes' : tuple(ModelPrefixes),
              'input_list': tuple(input_list),
              'multisolve' : multisolve
              }
        del SSAall
        RESULT_DATA.append(DD)
        RESULT_LPS.append(mGlp)


        FR = file(os.path.join(work_dir, '%s_rpt.txt' % MultiFNbase), 'w')
        FR.write('MultipleInput(%s states) Result:\n' % len(multisolve))
        #FR.write(' Flux diff: %s\n' % all_flux_diff.sum())
        FR.write(' optimum: %s = %s\n' % (mGobid, mGobjval))
        FR.write(' status: %s\n\n' % CBSolver.cplx_getSolutionStatus(mGlp))
        if options['algorithm'] == "ABS_L1_METHOD":
            FR.write('Distance: %s\n' % 'Sum(||J_a_x| - |J_b_x|| ...)')
        elif options['algorithm'] == "L1_METHOD":
            FR.write('Distance: %s\n' % 'sum(|J_a_x - J_b_x| ...)')
        FR.write('Big M: %s\n' % options['L1_BIGM_VALUE'])
        FR.write('With cost scaling: %s\n\n' % options['with_cost'])
        FR.write('MultipleInput(%s states) Solve Time\n' % len(multisolve))
        FR.write('Time: %s\nTime(m): %s\nTime(h): %s\n\n' % (TIME_END-TIME_START, (TIME_END-TIME_START)/60.0, (TIME_END-TIME_START)/60.0**2))

        FR.close()

        FD = file(os.path.join(work_dir, '%s_dat.test.csv' % MultiFNbase), 'w')
        head = 'Flux, '
        for p in ModelPrefixes:
            head += '%s,' % p
        FD.write('%sDiff,ProtCost,Diff(scaled)\n' % head)
        for F in range(FlxAO_scaled.shape[0]):
            FD.write('%s,' % FluxNamesOut[F])
            ##  for C in range(FlxAO_scaled.shape[1]):
                ##  if C != FlxAO_scaled.shape[1]-1:
                    ##  FD.write('%s,' % FlxAO_scaled[F, C])
                ##  else:
                    ##  FD.write('%s\n' % FlxAO_scaled[F, C])
            for C in range(FlxAO_scaled.shape[1]):
                FD.write('%s,' % FlxAO_scaled[F, C])
            ##  FD.write('%s,%s,%s\n' % (FluxScaling[F], FluxScaling[F]*FlxAO_scaled[F,-1]))
            FD.write('%s,%s\n' % (FluxScaling[F], FlxAO_scaled[F,-1]))
        FD.flush()
        FD.close()
        FD = file(os.path.join(work_dir, '%s_dat.csv' % MultiFNbase), 'w')
        head = 'Flux, '
        for p in ModelPrefixes:
            head += '%s,' % p
        FD.write('%sProtCost\n' % head)
        for F in range(FlxAO_J.shape[0]):
            FD.write('%s,' % FluxNamesOut[F])
            ##  for C in range(FlxAO_J.shape[1]):
                ##  if C != FlxAO_J.shape[1]-1:
                    ##  FD.write('%s,' % FlxAO_J[F, C])
                ##  else:
                    ##  FD.write('%s\n' % FlxAO_J[F, C])
            for C in range(FlxAO_J.shape[1]):
                FD.write('%s,' % FlxAO_J[F, C])
            FD.write('%s\n' % (FluxScaling[F]))
        FD.flush()
        FD.close()


    if options['with_cost']:
        OFname = '%s_P(%s)_C' % (options['input_model_file'], ModelPrefixStr)
    else:
        OFname = '%s_P(%s)' % (options['input_model_file'], ModelPrefixStr)
    if options['algorithm'] == "ABS_L1_METHOD":
        CBTools.storeObj(RESULT_DATA, os.path.join(work_dir, '%s_ABSL1_multi.dat' % OFname))
    elif options['algorithm'] == "L1_METHOD":
        CBTools.storeObj(RESULT_DATA, os.path.join(work_dir, '%s_L1_multi.dat' % OFname))

    ## multi file
    #if options['algorithm'] == "ABS_L1_METHOD":
        #FD = file(os.path.join(work_dir, '%s_ABSL1_multi.csv' % OFname), 'w')
    #elif options['algorithm'] == "L1_METHOD":
        #FD = file(os.path.join(work_dir, '%s_L1_multi.csv' % OFname), 'w')

    #dat_arr = [a['solution'] for a in RESULT_DATA]
    #dat_arr = numpy.hstack(dat_arr)

    #hcntr = 0
    #h2cntr = 0
    #head = 'Flux, '
    #for p in (ModelPrefixes)*len(dev_factor):
        #if hcntr <= len(ModelPrefixes)-1:
            #head += '%s(%s%s),' % (p, '', dev_factor[h2cntr])
            #hcntr += 1
        #else:
            #h2cntr += 1
            #head += '%s(%s%s),' % (p, '', dev_factor[h2cntr])
            #hcntr = 0

    #print '\n'
    #FD.write('%s\n' % head[:-1])
    #for F in range(dat_arr.shape[0]):
        #FD.write('%s,' % FluxNamesOut[F])
        #for C in range(dat_arr.shape[1]):
            #if C != dat_arr.shape[1]-1:
                #FD.write('%s,' % dat_arr[F, C])
            #else:
                #FD.write('%s\n' % dat_arr[F, C])
    #del dat_arr

    # start multi output
    dat_arr = numpy.hstack([a['solution'] for a in RESULT_DATA])
    if options['algorithm'] == "ABS_L1_METHOD":
        FDNAME = os.path.join(work_dir, '%s_ABSL1_multi.csv' % OFname)

    elif options['algorithm'] == "L1_METHOD":
        FDNAME = os.path.join(work_dir, '%s_L1_multi.csv' % OFname)

    H = []
    for R in range(len(dev_factor)):
        for S in range(len(ModelPrefixes)):
            H.append('%s(%s%s)' % (ModelPrefixes[S], '', dev_factor[R]))
    print(H)
    print(dat_arr.shape)
    CBTools.exportLabelledArrayWithHeader2CSV(dat_arr, fname=FDNAME, names=FluxNamesOut, header=H)
    del dat_arr
    # end multi output

    opt_out = [(a['dev_factor'],a['obj_val'],a['flux_diff']) for a in RESULT_DATA]
    if options['algorithm'] == "ABS_L1_METHOD":
        FR2 = file(os.path.join(work_dir, '%s_ABSL1_opt.csv' % OFname), 'w')
    elif options['algorithm'] == "L1_METHOD":
        FR2 = file(os.path.join(work_dir, '%s_L1_opt.csv' % OFname), 'w')
    head = 'ObjDist(%s),optimum (%s),raw diff (%s)' % (ModelPrefixStr,ModelPrefixStr,ModelPrefixStr)
    FR2.write('%s\n' % head)
    for d in opt_out:
        FR2.write('%s,%s,%s\n' % (d[0], d[1], d[2]))
    FR2.close()


    if dofva:
        MemDumpFiles = []
        for L in range(len(RESULT_LPS)):
            mdtempFile = os.path.join(work_dir,'_memDump_%s.lp' % L)
            RESULT_LPS[L].write(mdtempFile, filetype='lp')
            MemDumpFiles.append(mdtempFile)

        for L in range(len(RESULT_LPS)-1,-1,-1):
            tmp = RESULT_LPS.pop(L)
            del tmp
        del RESULT_LPS

        print(MemDumpFiles)

        fva_time_start = time.time()
        fva_meta_res = []

        if HAVE_XLWT:
            wbk = xlwt.Workbook()
            wbk2 = xlwt.Workbook()
            wbk3 = xlwt.Workbook()

        GO = True
        pCntr = [1 for a in range(len(RESULT_DATA[L]['prefixes']))]
        for L in range(len(MemDumpFiles)):
            if HAVE_XLWT:
                dfStr = str(RESULT_DATA[L]['dev_factor'])
                sheet = wbk.add_sheet(dfStr)
                sheet.write(0, 0, 'id(%s)' % dfStr)
                sheet.write(0, 1, 'J(%s)' % dfStr)
                sheet.write(0, 2, 'min(%s)' % dfStr)
                sheet.write(0, 3, 'max(%s)' % dfStr)

            dLP = CBSolver.cplex.Cplex(MemDumpFiles[L])
            # remove artificial variables
            lp_var = dLP.variables.get_names()
            if fva_selected_reactions == None:
                lp_var = [j for j in lp_var if 'xvar_' not in j]
                lp_var = [j for j in lp_var if 'absL_' not in j]
                lp_var = [j for j in lp_var if 'zvar1' not in j]
            else:
                for r in fva_selected_reactions:
                    assert r in lp_var, 'invalid input to FVA reaction selection'
                lp_var = fva_selected_reactions
            pFx = RESULT_DATA[L]['prefixes']
            print(pFx)
            if HAVE_XLWT:
                if GO:
                    sheetL3 = []
                    for p_ in pFx:
                        s0 = wbk3.add_sheet(p_)
                        s0.write(0, 0, p_)
                        s0.write(0, 1, 'MaxBio pc(%s)' % p_.replace('_',''))
                        s0.write(0, 2, 'opt(%s)' % p_.replace('_',''))
                        s0.write(0, 3, 'min(%s)' % p_.replace('_',''))
                        s0.write(0, 4, 'max(%s)' % p_.replace('_',''))
                    for v_ in lp_var:
                        s1 = wbk2.add_sheet(v_)
                        s1.write(0, 0, v_)
                        s1.write(0, 1, 'opt')
                        s1.write(0, 2, 'min')
                        s1.write(0, 3, 'max')
                    GO = False

                    #g_R_EX_etoh_e	Opt(g)	Min(g)	Max(g)

            FVAtol = None
            FVAoptper = 100.0
            if 'FVAtol' in options:
                FVAtol = options['FVAtol']
            if 'FVAoptper' in options:
                FVAoptper = options['FVAoptper']

            if options['procs'] == None:
                print('\nOld single threaded algorithm\n')
                fva_res, fva_names = CBSolver.cplx_MultiFluxVariabilityAnalysis(dLP, selected_reactions=lp_var, tol=FVAtol,\
                                                                                rhs_sense='upper', optPercentage=FVAoptper,\
                                                                                work_dir=work_dir, debug=False)
            else:
                print('\nNew multi-threaded algorithm using {} threads.\n'.format(options['procs']))
                fva_res, fva_names = CBMultiCore.runMultiCoreMultiEnvFVA(dLP, selected_reactions=lp_var, tol=FVAtol,\
                                                                         rhs_sense='upper', optPercentage=FVAoptper,\
                                                                         work_dir=work_dir, debug=False, procs=options['procs'])
            fva_meta_res.append((fva_res,fva_names))
            FVAOname = '%s(%s)_fva.csv' % (OFname,RESULT_DATA[L]['dev_factor'])
            FVAOname = os.path.join(work_dir, FVAOname)
            CBTools.exportLabelledArrayWithHeader(fva_res, fname=FVAOname, fva_names=fva_names, header=['Opt', 'Min', 'Max', 'MinS', 'MaxS'], sep=',', format='%f')
            if HAVE_XLWT:
                rfact = 8
                for j_ in range(1, len(fva_names)+1):
                    # ex 1
                    J = fva_names[j_-1]
                    sheet.write(j_, 0, J)
                    # ex 2
                    sheet2 = wbk2.get_sheet(j_-1)
                    sheet2.write(L+1, 0, RESULT_DATA[L]['dev_factor'])
                    sheet2.write(L+1, 1, fva_res[j_-1, 0])
                    res = fva_res[j_-1, 1]
                    if numpy.isnan(res):
                        #res = fva_res[j_-1, 0]
                        res = 'NaN'
                    else:
                        res = round(res, rfact)
                    sheet2.write(L+1, 2, res)
                    res = fva_res[j_-1, 2]
                    if numpy.isnan(res):
                        #res = fva_res[j_-1, 0]
                        res = 'NaN'
                    else:
                        res = round(res, rfact)
                    sheet2.write(L+1, 3, res)
                    for d_ in range(1,4):
                        res = fva_res[j_-1, d_-1]
                        if numpy.isnan(res):
                            #res = fva_res[j_-1, 0]
                            res = 'NaN'
                        else:
                            res = round(res, rfact)
                        sheet.write(j_, d_, res)
                    # ex 3

                    for p_ in range(len(pFx)):
                        if J.startswith(pFx[p_]):
                            sheet3 = wbk3.get_sheet(p_)
                            #wbk3.save('fuckme.xls')
                            sheet3.write(pCntr[p_], 0, J)
                            sheet3.write(pCntr[p_], 1, RESULT_DATA[L]['dev_factor'])
                            sheet3.write(pCntr[p_], 2, fva_res[j_-1, 0])
                            res = fva_res[j_-1, 1]
                            if numpy.isnan(res):
                                #res = fva_res[j_-1, 0]
                                res = 'NaN'
                            else:
                                res = round(res, rfact)
                            sheet3.write(pCntr[p_], 3, res)
                            res = fva_res[j_-1, 2]
                            if numpy.isnan(res):
                                #res = fva_res[j_-1, 0]
                                res = 'NaN'
                            else:
                                res = round(res, rfact)
                            sheet3.write(pCntr[p_], 4, res)
                            pCntr[p_] = pCntr[p_] + 1
            del dLP
        if HAVE_XLWT:
            t = str(time.time()).split('.')[0][-5:]
            xlname = 'data_%s_(%s)' % (OFname.replace('.xml','').replace('base',''), t)
            wbk.save(os.path.join(work_dir, xlname+'.fva_per_objdist.xls'))
            wbk2.save(os.path.join(work_dir, xlname+'.fva_per_flux.xls'))
            wbk3.save(os.path.join(work_dir, xlname+'.fva_per_input.xls'))
        fva_time_end = time.time()
        print('\n*****\nTotal ({} LP\'s) time for FVA ({} threads) min: {}\n*****\n'.format(len(MemDumpFiles)*len(lp_var), options['procs'], (fva_time_end-fva_time_start)/60.0))
        F = file(os.path.join(work_dir, '%s_fva_stats.txt' % OFname), 'a')
        F.write('*****\nTotal (%s LP\'s) time for FVA (%s threads) min: %s\n*****\n' % (len(MemDumpFiles)*len(lp_var), options['procs'], (fva_time_end-fva_time_start)/60.0))
        F.close()
        CBTools.storeObj(fva_meta_res, os.path.join(work_dir, '%s_fva_meta.dat' % OFname))
    return True
