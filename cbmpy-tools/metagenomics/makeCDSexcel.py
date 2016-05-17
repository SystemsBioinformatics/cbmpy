"""
Read in a full genbank file and extract all annotation. Returns three dictionaries containing CDS, Gene and Other information.

 - *fname* the full filename of the input filename in GenBank full format (typically *.gbk or *.gb)

"""
import Bio, os, time, json, xlwt
from Bio import SeqIO


fnames = ['NC_003030_1.gbk','NC_005090_1.gbk']

output = {}

for gbkf in fnames:
    output[gbkf] = []
    GBFile = file(gbkf, 'r')
    GBcds = Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
    cntr = 0
    #ltags = []
    for cds in GBcds:
        if cds.seq != None:
            #ltags.append(cds.name)
            data = {'pid' : cds.id,
                    'product' : cds.annotations['product'],
                    'seq' : str(cds.seq),
            }
        output[gbkf].append(data)
    GBFile.close()

wb = xlwt.Workbook(encoding='utf-8')
for f in output:
    ws = wb.add_sheet(f)
    for p in range(len(output[f])):
        ws.write(p, 0, output[f][p]['pid'])
        ws.write(p, 1, output[f][p]['product'])
        ws.write(p, 2, output[f][p]['seq'])

wb.save('{}.xls'.format('community_cds'))
