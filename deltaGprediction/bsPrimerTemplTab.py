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
argparser.add_argument('-t', type=str, default='bs_se',required=False, help='Original data of bam (bs_pe, bs_se, no_bs or rna)')
args = argparser.parse_args()

def make_templ_primer_dic(primerfile,templFasta, type):
    '''
    Makes dictionary of template-primer pairs. Takes template sequences defined in primedreg.fa and
    primer sequences is reverse complementary sequence of the first 6 bases of the untrimmed reads.
    Type is for matching the right read name between fasta files (bs_se adds comment, no_bs doesn't)
    '''
    templDic={}
    with ps.FastxFile(templFasta) as templ:
        for entry in templ:
            seq, name = entry.sequence.upper(), entry.name
            templDic[name]=[seq]
    with ps.FastxFile(primerfile) as f:
        for entry in f:
            if type=='bs_se':
                name = '_'.join([entry.name, entry.comment])
            elif type=='no_bs':
                name = entry.name
            else:
                print('SPECIFY TYPE!! bs_se or no_bs')
            if name in templDic.keys():
                templDic[name].append(entry.sequence[0:6])
    return(templDic)

fasta = args.primedreg
primerfile = args.primerinput
abundanceFile = args.abfile
type = args.t

templDic = make_templ_primer_dic(primerfile,fasta, type)

tabAb = pd.read_csv(abundanceFile, index_col=0, header=None)
df = make_occurrencies_tbl(templDic)
path = '/'.join(fasta.split('/')[:-1])
sample = fasta.split('/')[-1].split('.')[0]
df = fillNsortPTmatrix(df, tabAb)
df.to_csv(path + sample + '.ptCounts.qualFilt.parallel.csv')
