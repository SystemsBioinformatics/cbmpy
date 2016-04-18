import os, time, numpy, json, pprint, csv
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import cbmpy

iDir = os.path.join(cDir, 'input')
oDir = os.path.join(cDir, 'output')
tDir = os.path.join(cDir, 'temp')

for p in [oDir, tDir]:
    if not os.path.exists(p):
        os.makedirs(p)

reac_files = [f for f in os.listdir(iDir) if '(reactions)' in f]
metab_files = [f for f in os.listdir(iDir) if '(metabolites)' in f]
net_files = [f for f in os.listdir(iDir) if '(network_reactions)' in f]
gene_files = [f for f in os.listdir(iDir) if '(gene_assoc)' in f]

print(len(reac_files))
print(len(metab_files))
print(len(net_files))
print(len(gene_files))

assert len(reac_files) == len(metab_files) == len(net_files) == len(gene_files), '\nHave to have equal number of all files.'

reac_files.sort()
metab_files.sort()
net_files.sort()
gene_files.sort()

#print(reac_files)
#print(metab_files)
#print(net_files)
#print(gene_files)

#DEBUG
fmax = 1
reac_files = reac_files[0:fmax]
metab_files = metab_files[0:fmax]
net_files = net_files[0:fmax]
gene_files = gene_files[0:fmax]


print(reac_files)
print(metab_files)
print(net_files)
print(gene_files)
#DEBUG

for m_ in range(len(reac_files)):
    FR = open(os.path.join(iDir, reac_files[m_]), 'rU')
    FN = open(os.path.join(iDir, net_files[m_]), 'rU')
    FG = open(os.path.join(iDir, gene_files[m_]), 'rU')
    FM = open(os.path.join(iDir, metab_files[m_]), 'rU')

    FMtmp = open(os.path.join(iDir, 'temp_metab'), 'w')
    csvwrite = csv.writer(FMtmp, dialect='excel')
    for l in FM:
        L = l.split('\t')
        L = [i.strip() for i in L]
        csvwrite.writerow(L)
    FMtmp.flush()
    FMtmp.close()
    FM.close()
    del FMtmp, csvwrite
    FM = open(os.path.join(iDir, 'temp_metab'), 'r')

    Rtab = csv.reader(FR)
    Mtab = csv.reader(FM)
    Ntab = csv.reader(FN)
    Gtab = csv.reader(FG)

    MetD = {}
    Metk = None
    GO = True
    HEAD = True
    while GO:
        try:
            m = Mtab.next()
            if HEAD:
                Metk = m
                HEAD = False
            else:
                MetD[m[0]] = dict(zip(Metk, m))
                MetD[m[0]]['boundary'] = False
                MetD[m[0]]['compartment'] = 'cell'
        except StopIteration:
            GO = False

    ReacD = {}
    Reack = None
    GO = True
    HEAD = True
    rcntr = 0
    while GO:
        try:
            r = Rtab.next()
            if HEAD:
                Reack = r
                HEAD = False
            else:
                ReacD[rcntr] = dict(zip(Reack, r))
                ReacD[rcntr]['reversible'] = True
                ReacD[rcntr]['compartment'] = 'cell'
                ReacD[rcntr].pop('lowerbound')
                ReacD[rcntr].pop('upperbound')
                rcntr += 1
        except StopIteration:
            GO = False

    NetD = {}
    GO = True
    while GO:
        try:
            rid = Ntab.next()
            sub = Ntab.next()
            prod = Ntab.next()
            reagents = []
            for s in range(1, len(sub), 2):
                reagents.append((-int(sub[s]), sub[s+1]))
            for p in range(1, len(prod), 2):
                reagents.append((int(prod[p]), prod[p+1]))
            NetD[rid[0]] = reagents
        except StopIteration:
            GO = False

    for r in ReacD:
        if ReacD[r]['id'] in NetD:
            ReacD[r]['reagents'] = NetD[ReacD[r]['id']]

    GprD = {}
    GO = True
    while GO:
        try:
            l = Gtab.next()
            print(l)
            for t in range(len(l)):
                if t == 0:
                    key = l[t]
                    GprD[key] = '('
                elif l[t] == '':
                    pass
                elif l[t] == 'OR':
                    GprD[key] += ') or ('
                else:
                    glist = [ReacD[a]['NCBIaccession'] for a in ReacD if ReacD[a]['KEGG_ORTHOLOGY'] == l[t]]
                    print(glist)
                    if len(glist) == 0:
                        pass
                    elif len(glist) > 1:
                        gprm = '('+' or '.join(glist)+')'
                    else:
                        gprm = glist[0]
                    GprD[key] += '{} and '.format(gprm)
            GprD[key] += ')'
            GprD[key] = GprD[key].replace(' and )', ')')
        except StopIteration:
            GO = False

    # add raw OR gene associations to ReacD and overwrite it with predefined ones
    gass = {}
    for r in list(ReacD.keys()):
        if ReacD[r]['id'] not in gass:
            gass[ReacD[r]['id']] = {'gbank_seq_' + ReacD[r]['NCBIaccession'] : ReacD[r]['Sequence'],
                                    'organisms' : [ReacD[r]['Organism']]}
        else:
            gass[ReacD[r]['id']]['gbank_seq_' + ReacD[r]['NCBIaccession']] = ReacD[r]['Sequence']
            gass[ReacD[r]['id']]['organisms'].append(ReacD[r]['Organism'])
            ReacD.pop(r)
    for r in list(ReacD.keys()):
        newid = ReacD[r]['id']
        ReacD[newid] = ReacD.pop(r)
        ReacD[newid]['organisms'] = str(gass[newid].pop('organisms'))
        gpr = '('+' or '.join(gass[newid])+')'
        gpr = gpr.replace('gbank_seq_','')
        ReacD[newid]['GENE_ASSOCIATION'] = gpr
        for g in list(gass[newid]):
            ReacD[newid][g] = gass[newid].pop(g)

    for r in ReacD:
        if r in GprD:
            ReacD[r]['GENE_ASSOCIATION'] = GprD[ReacD[r]['id']]
        ReacD[r].pop('Sequence')
        ReacD[r].pop('NCBIaccession')
        ReacD[r].pop('Organism')
        ReacD[r].pop('Gene')

    F = open(os.path.join(tDir, reac_files[m_].replace('.csv','.reactions.json')),'w')
    json.dump(ReacD, F, indent=1, separators=(',', ': '))
    F.close()
    F = open(os.path.join(tDir, reac_files[m_].replace('.csv','.metabolites.json')),'w')
    json.dump(MetD, F, indent=1, separators=(',', ': '))
    F.close()

    print('\n\n**********\nCSV parsed into dictionaries, time to build a model\n**********\n\n')

    dmod = cbmpy.CBModelTools.quickDefaultBuild(model_name=reac_files[m_],\
                                         Reactions=ReacD,\
                                         Species=MetD,\
                                         Bounds={},
                                         Objective_function={},
                                         infinity=cbmpy.INF)

    for s_ in dmod.species:
        if s_.getAnnotation('chemformula') not in ['Not found', '']:
            s_.setChemFormula(s_.getAnnotation('chemformula'))
        s_.deleteAnnotation('chemformula')
        s_.deleteAnnotation('fixed')
        if s_.getAnnotation('charge') != '':
            s_.setCharge(float(s_.getAnnotation('charge')))
        s_.deleteAnnotation('charge')

    for r_ in dmod.reactions:
        r_.deleteAnnotation('compartment')


    cbmpy.writeSBML3FBC(dmod, reac_files[m_].replace('.csv','.xml'), directory=oDir,\
                        gpr_from_annot=False,
                        add_groups=False,
                        add_cbmpy_annot=True,
                        add_cobra_annot=False)
    cbmpy.writeModelToCOMBINEarchive(dmod, fname=reac_files[m_].replace('-(reactions).csv',''), directory=oDir,
                                    withExcel=True,
                                    vc_given='Brett',
                                    vc_family='Olivier',
                                    vc_email='b.g.olivier@vu.nl',
                                    vc_org='Vrije Universiteit Amsterdam',
                                    add_cbmpy_annot=True,
                                    add_cobra_annot=False)
    #del dmod, ReacD, NetD, MetD, GprD

