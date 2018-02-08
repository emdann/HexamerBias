from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd
import multiprocessing

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('bam', type=str, help='Bam input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
argparser.add_argument('ref', type=str, help='template or primer')
args = argparser.parse_args()

fasta = args.primedreg
bamfile=args.bam
ref = args.ref

# bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/zebrafish/gk2a-2.sam'
# fasta='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2_primed_seq.fa'

def makeTemplPrimerDic(bamfile,templFasta):
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name
            templDic[name]=[seq]
    with ps.AlignmentFile(bamfile,"rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.flag==0:
                templDic[r.qname].append(r.seq[0:6])
            elif r.flag==16:
                templDic.pop(r.qname)
    return(templDic)

def cellSpecificTbl(cellDic, cellsOI):
    '''
    Makes primer VS template occurrency table for a given subset of cells.
    Input: dict of {cell:{read:[template sequence, primer sequence]}}
    Output: dict of {cell name:matrix}
    '''
    tblCellDic={}
    for cell in CellsOI:
        tblCellDic[cell] = make_occurrencies_tbl(cellDic[cell])
    # ---> primers on columns, templates on rows <---
    return(tblCellDic)

templDic = makeTemplPrimerDic(bamfile,fasta)
cellDic={}
for name,seqs in templDic.items():
    cell = name.split(':')[-1]
    if cell not in cellDic.keys():
        cellDic[cell]={}
    cellDic[cell][name]=seqs
highcovCells = [i for i in cellDic.keys() if len(cellDic[i].values())>10000]
tblCellDic = cellSpecificTbl(templDic, highcovCells)

for cell in tblCellDic:
    tblCellDic[cell].to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/cell'+cell+'ptCounts.csv')
