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

def kmer_distTSS(params):
	# compute the distance from a given tss for each hexamer in the region flanking the tss
	# seq = sequence of the flanking region
	# 
	seq,k = params
	compl_seq=str(Seq(seq, generic_dna).complement())
	spl_seq=list(seq)
	# bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]	
	basesPlus=[[spl_seq[i],'T'] if ((spl_seq[i]=="C") & (spl_seq[i+1]=="G")) else spl_seq[i] for i in range(len(spl_seq)-1)]
	basesPlus=['T' if i=='C' else i for i in basesPlus]
	spl_seq_compl=list(compl_seq)
	basesMinus=[[spl_seq_compl[i],'T'] if ((spl_seq_compl[i]=="C") & (spl_seq[i-1]=="G")) else spl_seq_compl[i] for i in range(len(spl_seq_compl)-1)]
	basesMinus=['T' if i=='C' else i for i in basesMinus]
	tss=(len(seq)/2)-1
	tss_dist={}
	for i in range(0, len(basesPlus)-k+1):
	    hex = basesPlus[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if(''.join(hex_perm[x]) not in tss_dist.keys()):
	    		tss_dist[''.join(hex_perm[x])]=[]
	    	tss_dist[''.join(hex_perm[x])].append(i-tss)
	tss_distMinus={}
	for i in range(0, len(basesMinus)-k+1):
	    hex = basesMinus[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if(''.join(hex_perm[x]) not in tss_distMinus.keys()):
	    		tss_distMinus[''.join(hex_perm[x])]=[]
		    	tss_distMinus[''.join(hex_perm[x])].append(i-tss)	
	tss_distRevMinus={}
	for key,val in tss_distMinus.items():
		tss_distRevMinus[key[::-1]]=val
	for key,val in tss_distRevMinus.items():
		if key not in tss_dist.keys():
			tss_dist[key]=[]
		for pos in val:
			tss_dist[key].append(pos)
	return(tss_dist)

def make_occurrencies_tbl(tss_dist):
	dic_oc={}
	hexs=list(tss_dist.keys())
	counts=[collections.Counter(i) for i in tss_dist.values()]
	col_ind = [i for ids in tss_dist.values() for i in ids]
	col_ind = list(set(col_ind))
	for i in list(range(len(hexs))):
		dic_oc[hexs[i]]=counts[i]
	oc_tbl=pd.DataFrame(dic_oc).T
	oc_tbl=oc_tbl.fillna(0)
	return(oc_tbl)

# chromosome=fasta_file.split('/')[-1].split('.')[0]
chromosome=args.fasta.split('/')[-1].split('.')[0]
# print chr
# refgen = pd.read_csv(refgen_file, sep="\t", usecols=[0,2,3,4], header=0, dtype={4:int}) 
refgen = pd.read_csv(args.refgen, sep="\t", usecols=[0,2,3,4], header=0, dtype={4:int}) 
refgen=refgen[refgen.chrom==chromosome]
refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)
tss=refgen.txStart
# print refgen
# with ps.FastxFile(fasta_file) as chr:
with ps.FastxFile(args.fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

flank_wid=3000
seqs=[]
for i in (tss-1):
	start_pos=i-flank_wid
	end_pos=i+flank_wid
	small_seq=seq[start_pos:end_pos]
	seqs.append(small_seq)


workers = multiprocessing.Pool(10)
tss_dist={}
for dist in workers.imap_unordered(kmer_distTSS, [ (seqs[i],args.k) for i in list(range(len(seqs)))]):
	for key,val in dist.items():
		if key not in tss_dist.keys():
			tss_dist[key]=[]
		for pos in val:
			tss_dist[key].append(pos)

oc_tbl=make_occurrencies_tbl(tss_dist)
print("#hex\t"+'\t'.join([str(i) for i in list(oc_tbl)]))
for x in list(range(len(oc_tbl))):
	print(oc_tbl.index[x]+'\t'+'\t'.join([str(i) for i in oc_tbl.ix[x]]))
















