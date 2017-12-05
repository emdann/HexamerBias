import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
import multiprocessing
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count hexamers in BS converted fasta file.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('cov2c', type=str, help='cytosine report input')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

def weightCountKmers(params):
	seq,df,k = params
	spl_seq=list(seq)
	probs=[[df.frac[df.pos==i+1].values[0], 1-df.frac[df.pos==i+1].values[0]] if i in list(df["pos"]-1) else [1] for i in list(range(len(spl_seq)))]
	# print("Probs computed.")
	## strand specific
	bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]
	# print("Bases computed.")
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

def strandSpecificCount(params):
	seq,df,k = params
	dfPlus = df[df.strand=="+"]
	dfMinus = df[df.strand=="-"]
	my_dna = Seq(seq, generic_dna)
	kmerCounts=collections.Counter()
	kmerCountsPlus=weightCountKmers((seq,dfPlus,k))
	kmerCountsMinus=weightCountKmers((seq,dfMinus,k))
	kmerCountsRevMinus=collections.Counter()
	for key,val in kmerCountsMinus.items():
		kmerCountsRevMinus[str(Seq(key, generic_dna).reverse_complement())]+=val
	kmerCounts+=kmerCountsPlus
	kmerCounts+=kmerCountsRevMinus
	return(kmerCounts)

# Load cov2c
# cov2c file processed with added fraction of methylation
# zcat cov2c_smp.deduplicated.bismark.cov.gz | awk '$4+$5!=0{print $N"\t"$4/($4+$5)}'
cov2c = pd.read_csv(args.cov2c, sep="\t", header=None) 
cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
cov2c=cov2c[-np.isnan(cov2c.frac)]
# print("cov2c loaded!")

# Load fasta file
with ps.FastxFile(args.fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()
# print("Fasta loaded!")

covs=[]
seqs=[]

bin_len=1000
start_pos=0
end_pos=bin_len
track=1
while end_pos<len(seq):
	# print("Binning bin no."+str(track))
	small_seq=seq[start_pos:end_pos]
	small_cov = cov2c[(cov2c.pos<end_pos) & (cov2c.pos>start_pos)]
	small_cov = small_cov.assign(pos=small_cov.pos-start_pos)
	seqs.append(small_seq)
	covs.append(small_cov)
	# weightCountKmers(small_seq,small_cov,6)
	start_pos = end_pos-5
	end_pos = start_pos + bin_len
	# track+=1

# print("Binning last bin!")
last_seq = seq[start_pos:]
last_cov = cov2c[(cov2c.pos < end_pos) & (cov2c.pos > start_pos)]
seqs.append(last_seq)
covs.append(last_cov)

workers = multiprocessing.Pool(10)
counts=collections.Counter()

# track=1
for kmerCounts in workers.imap_unordered(strandSpecificCount, [ (seqs[i],covs[i],args.k) for i in list(range(len(seqs)))]):
        # print("Adding count no. "+str(track))
        counts+=kmerCounts
        # track+=1

for kmer, abundance in counts.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")

