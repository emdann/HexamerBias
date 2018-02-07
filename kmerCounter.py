#!/usr/bin/env python3
import pysam
import argparse
import collections
import multiprocessing
import pandas as pd
# argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count kmers in fasta file. By Buys de Barbanson")
# argparser.add_argument('fasta', type=str, help='Fasta input')
# argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
# argparser.add_argument('-t', type=int, default=8, required=False, help='Amount of threads to use.')
# args = argparser.parse_args()

def find_kmers(params):
    string, k, addTo = params
    kmerCounts = collections.Counter() if addTo is None else addTo
    n = len(string)
    for i in range(0, n-k+1):
        kmerCounts[ string[i:i+k] ] += 1
    return kmerCounts

def cellKmersAbundance(params):
    '''
    Take dictionary of abundances counts per entry of fasta file and cell transcript count table and cell name
    return count based on number of transcript for that cell
    '''
    CountDic,countTranscr,cellname=params
    cell=countTranscr[cellname]
    umiCountDic=collections.Counter()
    for gene,counter in CountDic.items():
        if not cell[cell.index==gene].empty:
            n = int(cell[cell.index==gene])
            # print('gene ' + gene + ' is present ' + str(n) + ' times...')
            geneNewDict = collections.Counter()
            for hex,count in counter.items():
                geneNewDict[hex]+=count*n
            umiCountDic+=geneNewDict
    return({cellname:umiCountDic})


coutt='/hpc/hub_oudenaarden/aalemany/emma-adi/zebrafish/gk2a-2.coutt.csv'
fasta = '/hpc/hub_oudenaarden/abarve/genomes/Danio_rerio_Zv9_ens74_extended3_genes_ERCC92_GFPmod_geneids.fa'

## Read files
countT = pd.read_csv(coutt, sep='\t', index_col=0)
countDic={}
with ps.FastxFile(fasta) as f:
    for entry in f:
        countDic[entry.name]=find_kmers((entry.sequence.upper(), 6, None))

# cellDic={}
# for cell in countT:
#     print('cell no. '+str(cell))
#     cellDic[cell]=cellKmersAbundance((countDic,countT[cell]))
#     print('done')

## Run on multiple cores
workers = multiprocessing.Pool(8)
finalKmerCounts = collections.Counter()

cellDic={}
for cell,counter in workers.imap_unordered(cellKmersAbundance, [ (countDic, countT, cell) for cell in countT]):
    cellDic[cell]=counter

outputTab = pd.DataFrame.from_dict(cellDic).to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2.cellAbundance.csv')
