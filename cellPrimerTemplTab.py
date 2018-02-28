from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd
import multiprocessing
import os

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts tables (run on 10 cores) \n By Emma Dann")
argparser.add_argument('bam', type=str, help='Bam input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
argparser.add_argument('-t', type=str, default='rna', required=False, help='Type of bamfile')
args = argparser.parse_args()

def make_templ_primer_dic(bamfile,templFasta, type='rna'):
    '''
    Makes dictionary of template-primer pairs, removing pairs if the phred score is low
    '''
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name
            templDic[name]=[seq]
    with ps.AlignmentFile(bamfile,"rb") as bam:
        if type=='rna':
            for r in bam.fetch(until_eof=True):
                if r.flag==0 and not any(q<32 for q in r.query_qualities[:6]):
                    templDic[r.qname].append(r.seq[0:6])
                elif r.flag==0 and any(q<32 for q in r.query_qualities[:6]):
                    templDic.pop(r.qname)
        if type=='bs':
            for r in bam.fetch(until_eof=True):
                if r.is_read1:
                    templDic[r.qname].append(r.seq[0:6])
    return(templDic)

def make_cell_pt_table(params):
    '''
    Makes primer VS template occurrency table, filling missing values and
    sorting in alphabetical order.
    Input: dict of {read:[template sequence, primer sequence]}
    Output: matrix
    '''
    ptDic,cellAb,cellname = params
    cellMatrix = fillNsortPTmatrix(make_occurrencies_tbl(ptDic), cellAb)
    # ---> primers on columns, templates on rows <---
    return({cellname:cellMatrix})

def split_pt_dic(templDic):
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
    return(cellDic)

def save_ptCounts(cellDic,cellsOI,fasta, cores=10):
    abundanceFile = fasta.strip('.primedreg.fa')+'.cellAbundance.csv'
    tabAb = pd.read_csv(abundanceFile, index_col=0)
    path = '/'.join(fasta.split('/')[:-1])
    if path:
        outpath = path + '/ptCounts/'
    else:
        outpath = './ptCounts/'
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    sample = bamfile.split('/')[-1].split('.')[0]
    workers = multiprocessing.Pool(cores)
    for cellMatrix in workers.imap_unordered(make_cell_pt_table, [ (cellDic[cell], tabAb[cell], cell) for cell in cellsOI]):
        for cell,matrix in cellMatrix.items():
            matrix.to_csv(outpath + sample + '.cell' + cell + '.ptCounts.qualFilt.parallel.csv')
    return('')

fasta = args.primedreg
bamfile = args.bam
type = args.type

templDic = make_templ_primer_dic(bamfile,fasta, type=type)
if type=='rna':
    cellDic = split_pt_dic(templDic)
    highcovCells = [i for i in cellDic.keys() if len(cellDic[i].values()) > 10000]
    save_ptCounts(cellDic, highcovCells, fasta)

if type=='bs':
    # abundanceFile = fasta.strip('.primedreg.fa')+'.cellAbundance.noN.csv'
    abundanceFile='mm10.cellAbundance.noN.csv'
    tabAb = pd.read_csv(abundanceFile, index_col=0, header=None)
    df = make_occurrencies_tbl(templDic)
    path = '/'.join(fasta.split('/')[:-1])
    sample = bamfile.split('/')[-1].split('.')[0]
    df = fillNsortPTmatrix(df, tabAb)
    df.to_csv(path + sample + '.ptCounts.qualFilt.parallel.csv')
