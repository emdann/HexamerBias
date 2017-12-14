import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
import scipy.sparse as sp
import multiprocessing
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
# from freeEnergy import compute_deltaG


argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='fasta input')
# argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
argparser.add_argument('bed', type=str, help='RefGen file of chr of interest (downloaded from UCSC genome browser)')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()


def compute_deltaG(hex):
	deltaG_tbl=np.array([-1.28,-1.54,-1.12,-1.72,-1.58,-2.07,-1.72,-2.53, -0.85, -1.73,-1.28,-1.58, -1.74,-2.49,-1.54,-2.07]).reshape((4,4))
	#hex="TTCGAT"
	pos={"A":0,"G":1,"T":2,"C":3}
	hex_G=0
	for i in list(range(len(hex)-1)):
		#print(deltaG_tbl[pos[hex[i]], pos[hex[i+1]]])
		hex_G += deltaG_tbl[pos[hex[i]], pos[hex[i+1]]]
	# hex_G += deltaG_tbl[pos[hex[-1]],].mean()
	return(hex_G)


def convert_seq(seq):
	spl_seq=list(seq)
	bases=[[spl_seq[i],'T'] if ((spl_seq[i]=="C") & (spl_seq[i+1]=="G")) else spl_seq[i] for i in range(len(spl_seq)-1)]
	bases.append(spl_seq[-1])
	bases=['T' if i=='C' else i for i in bases]
	return(bases)

def deltaG_pos(bases,k):
	tss=(len(bases)/2)-1
	pos_dic={}
	for i in range(0, len(bases)-k+1):
	    hex = bases[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if (i-tss) not in pos_dic.keys():
	    		pos_dic[i-tss]=[]
	    	# print ''.join(hex_perm[x])
	    	pos_dic[i-tss].append(compute_deltaG(''.join(hex_perm[x])))
	return(pos_dic)

def dG_distTSS(params, conv=True):
	# compute the distance from a given tss for each hexamer in the region flanking the tss
	# seq = sequence of the flanking region
	# 
	seq,k = params
	plStrand=seq
	minStrand=str(Seq(seq, generic_dna).reverse_complement())
	if conv:
		basesPlus = convert_seq(plStrand)
		basesMinus = convert_seq(minStrand)
	else:
		basesPlus=list(plStrand)
		basesMinus=list(minStrand)
	dg_distPlus = deltaG_pos(basesPlus,k)
	dg_distMinus = deltaG_pos(basesMinus,k)	
	for key,val in dg_distMinus.items():
		if key not in dg_distPlus.keys():
			dg_distPlus[key]=[]
		for pos in val:
			dg_distPlus[key].append(pos)
	return(dg_distPlus)


fasta=args.fasta
bedfile=args.bed

# fasta="/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa"
# bedfile="/hpc/hub_oudenaarden/edann/hexamers/rand_tss.mm10.txt"

bed = pd.read_csv(bedfile, sep="\t", header=None, usecols=[0,1,2], dtype={4:int}, names=["chrom", "start", "end"])
chrs = bed.chrom.unique()
genome={}
with ps.FastxFile(fasta) as gen:
	for chr in gen:
		if chr.name in chrs:
			genome[chr.name]=chr.sequence.upper()

flank_wid=3000
seqs=[]
for chr in chrs:
	refgen = bed[bed.chrom==chr]
	# print(refgen)
	# refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)
	tss = refgen.start
	tss = tss.drop_duplicates(keep='first', inplace=False)
	seq=genome[chr]
	for i in (tss-1):
		if i>flank_wid and (i+flank_wid)<len(seq):
			start_pos=i-flank_wid
			end_pos=i+flank_wid
		else:
			next
		small_seq=seq[start_pos:end_pos]
		seqs.append(small_seq)

dg= {}
workers = multiprocessing.Pool(8)
for dist in workers.imap_unordered(dG_distTSS, [ (seqs[i],6) for i in list(range(len(seqs)))]):
	for key,val in dist.items():
		if key not in dg.keys():
			dg[key]=[]
		dg[key].extend(val)

dgMeans={}
for k,v in dg.items():
	dgMeans[k]=np.mean(v)

for pos, meandg in dgMeans.items(): 
	print(f"{pos}\t{meandg}")

