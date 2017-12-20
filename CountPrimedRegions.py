import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
import multiprocessing
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count hexamers in BS converted fasta file.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('cov2c', type=str, help='cytosine report input')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

bedfile='/hpc/hub_oudenaarden/edann/hexamers/L1_primed_reg.bed'
fasta=args.fasta
cov2cfile=args.cov2c

def weightCountKmers(params):
	seq,df,k = params
	spl_seq=[['T','C'] if i == "C" else i for i in list(seq)]
	probs=[[1-df.frac[df.pos==i+1].values[0], df.frac[df.pos==i+1].values[0]] if i in list(df["pos"]-1) else [0.99,0.01] if len(spl_seq[i])==2 else [1] for i in list(range(len(spl_seq)))]
	# bases=[[spl_seq[i],'C'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'G'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]
	bases=spl_seq
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
	kmerCountsMinus=weightCountKmers((str(Seq(seq, generic_dna).complement()),dfMinus,k))
	kmerCountsRevMinus=collections.Counter()
	for key,val in kmerCountsMinus.items():
		kmerCountsRevMinus[key[::-1]]+=val
	kmerCounts+=kmerCountsPlus
	kmerCounts+=kmerCountsRevMinus
	return(kmerCounts)

def count_amplification(kmerCounts):
	countsAmp=collections.Counter()
	for hex,count in kmerCounts.items():
		countsAmp[hex]+=count
		countsAmp[str(Seq(hex, generic_dna).reverse_complement())]+= count
	return(countsAmp)


cov2c = pd.read_csv(cov2cfile, sep="\t", header=None) 
cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
cov2c=cov2c[-np.isnan(cov2c.frac)]
# print("cov2c loaded!")

chromosome = fasta.split('/')[-1].split('.')[0]
bed = pd.read_csv(bedfile, sep="\t", header=None, dtype={4:int}, names=["chrom", "start", "end"])
refgen = bed[bed.chrom==chromosome]

with ps.FastxFile(fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

# hexs=refgen.hex
covs=[]
seqs=[]
start = refgen.start
for i in start:
	start_pos = i
	end_pos = i+6
	small_cov = cov2c[(cov2c.pos<end_pos) & (cov2c.pos>start_pos)]
	small_cov = small_cov.assign(pos=small_cov.pos-start_pos)
	small_seq=seq[start_pos:end_pos]
	seqs.append(small_seq)
	covs.append(small_cov)


workers = multiprocessing.Pool(10)
counts=collections.Counter()

# track=1
for kmerCounts in workers.imap_unordered(strandSpecificCount, [ (seqs[i],covs[i],args.k) for i in list(range(len(seqs)))]):
        # print("Adding count no. "+str(track))
        counts+=kmerCounts
        # track+=1

for kmer, abundance in counts.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")

