from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd
import multiprocessing
import os
from Bio.Seq import Seq,MutableSeq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts tables for bs-seq \n By Emma Dann")
argparser.add_argument('primerinput', type=str, help='Untrimmed fasta input')
argparser.add_argument('primedreg', type=str, help='Fasta input of primed region sequences')
argparser.add_argument('abfile', type=str, help='Csv file of kmer abundance on reference genome')
args = argparser.parse_args()

def make_templ_primer_dic(primerfile,templFasta):
    '''
    Makes dictionary of template-primer pairs. Takes template sequences defined in primedreg.fa and
    primer sequences is reverse complementary sequence of the first 6 bases of the untrimmed reads.
    '''
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name
            templDic[name]=[seq]
    with ps.FastxFile(primerfile) as f:
        for entry in f:
            nameWbarcode = '_'.join([entry.name, entry.comment])
            if nameWbarcode in templDic.keys():
                templDic[nameWbarcode].append(entry.sequence[0:6])
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
