from __future__ import division

import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
import scipy.sparse as sp
import pickle


argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute avg base composition for n sequences of the same length from fasta file By Emma Dann")
argparser.add_argument('fasta', type=str, help='fasta input')
args = argparser.parse_args()

fasta=args.fasta
# fasta="/home/emma/tss_danRer.fa"

seq_dic={}
with ps.FastxFile(fasta) as chr:
 	for entry in chr:
 		seq_dic[entry.name]=entry.sequence.upper()

nuc_dic={"A":[], "C":[], "G":[], "T":[]}
for pos in range(len(seq_dic.values()[0])):
	nuc_dic["A"].append([i[pos] for i in seq_dic.values()].count("A")/len(seq_dic))
	nuc_dic["T"].append([i[pos] for i in seq_dic.values()].count("T")/len(seq_dic))
	nuc_dic["G"].append([i[pos] for i in seq_dic.values()].count("G")/len(seq_dic))
	nuc_dic["C"].append([i[pos] for i in seq_dic.values()].count("C")/len(seq_dic))

pd.DataFrame(nuc_dic).to_csv("nuc_"+fasta.split("/")[-1].strip(".fa")+".csv")