from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd
import multiprocessing
import os

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts tables (run on 10 cores) \n By Emma Dann")
argparser.add_argument('primerinput', type=str, help='Untrimmed fasta input')
argparser.add_argument('primedreg', type=str, help='Fasta input of primed region sequences')
argparser.add_argument('abfile', type=str, help='Csv file of kmer abundance on reference genome')
args = argparser.parse_args()

def make_templ_primer_dic(primerfile,templFasta, type='rna'):
    '''
    Makes dictionary of template-primer pairs, removing pairs if the phred score is low
    '''
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name.split('_')[0]
            templDic[name]=[seq]
    with ps.FastxFile(primerfile) as f:
        for entry in f:
            if entry.name in templDic.keys():
                templDic[entry.name].append(entry.sequence[0:6])
    return(templDic)

fasta = args.primedreg
primerfile = args.primerinput
abundanceFile = args.abfile

templDic = make_templ_primer_dic(primerfile,fasta)

tabAb = pd.read_csv(abundanceFile, index_col=0, header=None)
df = make_occurrencies_tbl(templDic)
path = '/'.join(fasta.split('/')[:-1])
sample = fasta.split('/')[-1].split('.')[0]
df = fillNsortPTmatrix(df, tabAb)
df.to_csv(path + sample + '.ptCounts.qualFilt.parallel.csv')
