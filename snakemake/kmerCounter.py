import pysam as ps
import argparse
import collections
import multiprocessing
import pandas as pd

def find_kmers(params):
    '''
    Count kmer abundance in sequence
    '''
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
