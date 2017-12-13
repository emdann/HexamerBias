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
import random
import multiprocessing
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='chr.fasta input')
# argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
argparser.add_argument('refgen', type=str, help='RefGen file of chr of interest (downloaded from UCSC genome browser)')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()


def convert_seq(seq):
	spl_seq=list(seq)
	bases=[[spl_seq[i],'T'] if ((spl_seq[i]=="C") & (spl_seq[i+1]=="G")) else spl_seq[i] for i in range(len(spl_seq)-1)]
	bases.append(spl_seq[-1])
	bases=['T' if i=='C' else i for i in bases]
	return(bases)

def kmer_pos(bases,k):
	tss=(len(bases)/2)-1
	tss_dist={}
	for i in range(0, len(bases)-k+1):
	    hex = bases[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if(''.join(hex_perm[x]) not in tss_dist.keys()):
	    		tss_dist[''.join(hex_perm[x])]=[]
	    	# print ''.join(hex_perm[x])
	    	tss_dist[''.join(hex_perm[x])].append(i-tss)
	return(tss_dist)

def kmer_distTSS(params, conv=False):
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
	tss_distPlus = kmer_pos(basesPlus,k)
	tss_distMinus = kmer_pos(basesMinus,k)	
	for key,val in tss_distMinus.items():
		if key not in tss_distPlus.keys():
			tss_distPlus[key]=[]
		for pos in val:
			tss_distPlus[key].append(pos)
	return(tss_distPlus)

def make_occurrencies_tbl(tss_dist):
	dic_oc={}
	hexs=list(tss_dist.keys())
	counts=[collections.Counter(i) for i in tss_dist.values()]
	for i in list(range(len(hexs))):
		dic_oc[hexs[i]]=counts[i]
	oc_tbl=pd.DataFrame(dic_oc).T
	oc_tbl=oc_tbl.fillna(0)
	return(oc_tbl)

chromosome=args.fasta.split('/')[-1].split('.')[0]
refgen = pd.read_csv(args.refgen, sep="\t", usecols=[0,2,3,4], header=0, dtype={4:int}) 
refgen = refgen[refgen.chrom==chromosome]
refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)
tss = refgen.txStart

# refgen=pd.read_csv(args.refgen, sep="\t", header=False, dtype={1:int,2:int}, names=["chrom", "start", "end"])
# refgen = refgen[refgen.chrom==chromosome]

with ps.FastxFile(args.fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

# seqs=[]
# for i in range(len(refgen)):
# 	start_pos=list(refgen.start)[i]
# 	end_pos=list(refgen.end)[i]
# 	small_seq=seq[start_pos:end_pos]
# 	seqs.append(small_seq)

flank_wid=3000
seqs=[]
for i in (tss-1):
	start_pos=i-flank_wid
	end_pos=i+flank_wid
	small_seq=seq[start_pos:end_pos]
	seqs.append(small_seq)


workers = multiprocessing.Pool(8)
tss_dist={}
for dist in workers.imap_unordered(kmer_distTSS, [ (seqs[i],args.k) for i in list(range(len(seqs)))]):
	for key,val in dist.items():
		if key not in tss_dist.keys():
			tss_dist[key]=[]
		for pos in val:
			tss_dist[key].append(pos)

oc_tbl=make_occurrencies_tbl(tss_dist)
print("hex\t"+'\t'.join([str(i) for i in list(oc_tbl)]))
for x in list(range(len(oc_tbl))):
	print(oc_tbl.index[x]+'\t'+'\t'.join([str(i) for i in oc_tbl.ix[x]]))
















