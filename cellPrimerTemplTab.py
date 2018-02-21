from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd
import multiprocessing

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts table \n By Emma Dann")
argparser.add_argument('bam', type=str, help='Bam input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
args = argparser.parse_args()

fasta = args.primedreg
bamfile = args.bam

def make_templ_primer_dic(bamfile,templFasta):
    '''
    Makes dictionary of template-primer pairs, removing pairs if the phred score is low
    '''
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name
            templDic[name]=[seq]
    with ps.AlignmentFile(bamfile,"rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.flag==0 and not any(q<32 for q in r.query_qualities[:6]):
                templDic[r.qname].append(r.seq[0:6])
            elif r.flag==0 and any(q<32 for q in r.query_qualities[:6]):
                templDic.pop(r.qname)
            # elif r.flag==16:
            #     templDic.pop(r.qname)
    return(templDic)

def cellSpecificTbl(cellDic, cellsOI):
    '''
    Makes primer VS template occurrency table for a given subset of cells.
    Fills missing values and sorts in alphabetical order
    Input: dict of {cell:{read:[template sequence, primer sequence]}}
    Output: dict of {cell name:matrix}
    '''
    tblCellDic={}
    for cell in cellsOI:
        tblCellDic[cell] = fillNsortPTmatrix(make_occurrencies_tbl(cellDic[cell]))
    # ---> primers on columns, templates on rows <---
    return(tblCellDic)


def make_cell_tp_dic(templDic):
    '''
    Split dictionary of reads:[template,primer] by cell.
    Output: dictionary of pt tables, one for each cell with more than 10k aligned reads
    '''
    cellDic={}
    for name,seqs in templDic.items():
        cell = name.split(':')[-1]
        if cell not in cellDic.keys():
            cellDic[cell]={}
        cellDic[cell][name]=seqs
    highcovCells = [i for i in cellDic.keys() if len(cellDic[i].values()) > 10000]
    tblCellDic = cellSpecificTbl(cellDic, highcovCells)
    return(tblCellDic)

templDic = make_templ_primer_dic(bamfile,fasta)
tblCellDic = make_cell_tp_dic(templDic)

outpath = '/'.join(fasta.split('/')[:-1]) + '/ptCounts/'
if not os.path.exists(outpath):
    os.makedirs(outpath)
sample = bamfile.split('/')[-1].split('.')[0]

for cell in tblCellDic:
    countsMatrix.to_csv(outpath + sample + '.cell' + cell + '.ptCounts.qualFilt.csv')
