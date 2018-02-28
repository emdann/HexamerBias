#!/usr/bin/env python3
import pysam as ps
import argparse
import collections
import multiprocessing
import pandas as pd
argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count kmers in fasta file. By Buys de Barbanson")
argparser.add_argument('coutc', type=str, help='UMI counts')
argparser.add_argument('refgen', type=str, help='Fasta file of reference genome')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
argparser.add_argument('-t', type=int, default=8, required=False, help='Amount of threads to use.')
argparser.add_argument('-o', type=str, required=False, help='Output directory')
args = argparser.parse_args()

def find_kmers(params):
    string, k, addTo = params
    kmerCounts = collections.Counter() if addTo is None else addTo
    n = len(string)
    for i in range(0, n-k+1):
        kmerCounts[ string[i:i+k] ] += 1
    return kmerCounts

def count_kmers_refgen(refgenFasta):
    countDic = {}
    with ps.FastxFile(refgenFasta) as f:
        for entry in f:
            countDic[entry.name]=find_kmers((entry.sequence.upper(), 6, None))
    return(countDic)

def cellKmersAbundance(params):
    '''
    Take dictionary of abundances counts per entry of fasta file and cell transcript count table and cell name
    return count based on number of transcript for that cell
    '''
    CountDic,countTranscr,cellname=params
    cell = countTranscr[cellname]
    umiCountDic = collections.Counter()
    for gene,counter in CountDic.items():
        if not cell[cell.index==gene].empty:
            n = int(cell[cell.index==gene])
            # print('gene ' + gene + ' is present ' + str(n) + ' times...')
            geneNewDict = collections.Counter()
            for hex,count in counter.items():
                geneNewDict[hex]+=count*n
            umiCountDic+=geneNewDict
    return({cellname:umiCountDic})

coutc = args.coutc
fasta = args.refgen

## Read files
countT = pd.read_csv(coutc, sep='\t', index_col=0)
countDic = count_kmers_refgen(fasta)

## Run on multiple cores
workers = multiprocessing.Pool(8)
finalKmerCounts = collections.Counter()

cellDic={}
for cellCounter in workers.imap_unordered(cellKmersAbundance, [ (countDic, countT, cell) for cell in countT]):
    if len(cellCounter.keys())==1:
        for cell,counter in cellCounter.items():
            cellDic[cell]=counter
    else:
        print('Something wrong...')

sample = coutc.split('/')[-1].split('.coutc')[0]
if args.o:
    outpath = args.o
else:
    outpath = '/'.join(fasta.split('/')[:-1]) + '/'

ab = pd.DataFrame.from_dict(cellDic)
noN = ab.T[[i for i in ab.index if 'N' not in i and 'Y' not in i]].T
print(noN)
outputTab = noN.to_csv(outpath + '/' + sample +'.cellAbundance.noN.csv')
