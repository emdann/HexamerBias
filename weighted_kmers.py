import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
# print("packages loaded!")

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count hexamers in BS converted fasta file.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('cov2c', type=str, help='cytosine report input')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

def weightCountKmers(seq,df,k):
	spl_seq=list(seq)
	probs=[[df.frac[df.pos==i+1].values[0], 1-df.frac[df.pos==i+1].values[0]] if i in list(df["pos"]-1) else [1] for i in list(range(len(spl_seq)))]
	## strand specific
	bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]
	kmerCounts = collections.Counter() 
	ls = []
	for i in range(0, len(bases)-k+1):
	    hex = bases[i:i+k]
	    hex_prob = probs[i:i+k]
	    hex_perm = list(it.product(*hex))
	    hex_prob_perm = list(it.product(*hex_prob))
	    # print(f"Processing window {i}")
	    for i in range(len(hex_perm)):
	    	kmerCounts[''.join(hex_perm[i])] += np.prod(np.array(hex_prob_perm)[i])
	return(kmerCounts)

# Load cov2c
# cov2c file processed with added fraction of methylation
# zcat cov2c_smp.deduplicated.bismark.cov.gz | awk '$4+$5!=0{print $N"\t"$4/($4+$5)}'
cov2c = pd.read_csv(args.cov2c, sep="\t", header=None) 
cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank", "frac"]

# Load fasta file
with ps.FastxFile(args.fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

counts=weightCountKmers(seq,cov2c,args.k)

for kmer, abundance in counts.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")






