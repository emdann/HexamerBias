import pysam as ps
import pandas as pd


bamfile='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2.bam'

geneid = 'ENSDARG00000009087_cd74a'


def makeBedForAllGenes(geneid, bamfile):
geneDic={}
with ps.AlignmentFile(bamfile,"rb") as bam:
    for r in bam.fetch(until_eof=True):
        if r.flag==0:
            gene = r.rname
            cell = r.qname.split(':')[-1]
            if gene not in geneDic.keys():
                geneDic[gene]={}
            if cell not in geneDic[gene].keys():
                geneDic[gene][cell]=[]
            geneDic[gene][cell].append([r.pos, r.pos+r.query_alignment_length])

def makeBedForGeneID(geneid, bamfile):
geneDic={}
with ps.AlignmentFile(bamfile,"rb") as bam:
    for r in bam.fetch(until_eof=True):
        if r.flag==0 and r.rname==geneid:
            gene = r.rname
            cell = r.qname.split(':')[-1]
            if cell not in geneDic.keys():
                geneDic[cell]=[]
            geneDic[cell].append([r.pos, r.pos+r.query_alignment_length])

with open(outputfile, 'wb') as out:
    print('cell', 'start', 'end', sep='\t', file=out)
    for k,v in geneDic.items():
        for val in v:
            print(k,val[0], val[1], sep='\t', file=out)
