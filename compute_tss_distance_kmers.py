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

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='chr.fasta input')
# argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
argparser.add_argument('refgen', type=str, help='RefGen file of chr of interest (downloaded from UCSC genome browser)')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

def kmer_distTSS(params):
	# compute the distance from a given tss for each hexamer in the region flanking the tss
	# seq = sequence of the flanking region
	# df = cov2c table of methylation fraction for the flanking region
	# 
	seq,k = params
	spl_seq=list(seq)
	# bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]	
	bases=[[spl_seq[i],'T'] if ((spl_seq[i]=="C") & (spl_seq[i+1]=="G")) else [spl_seq[i],'A'] if ((spl_seq[i]=="G") & (spl_seq[i-1]=="C")) else spl_seq[i] for i in range(len(spl_seq)-1)]
	# print(bases)
	tss=(len(seq)/2)-1
	tss_dist={}
	for i in range(0, len(bases)-k+1):
	    hex = bases[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if(''.join(hex_perm[x]) not in tss_dist.keys()):
	    		tss_dist[''.join(hex_perm[x])]=[]
	    	tss_dist[''.join(hex_perm[x])].append(i-tss)
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

# # cov2c = pd.read_csv(cov2c_file, sep="\t", header=None) 
# cov2c = pd.read_csv(args.cov2c, sep="\t", header=None) 
# cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
# cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
# cov2c=cov2c[-np.isnan(cov2c.frac)]

flank_wid=3000
# covs=[]
seqs=[]
for i in (tss-1):
	start_pos=i-flank_wid
	end_pos=i+flank_wid
	small_seq=seq[start_pos:end_pos]
	# small_cov = cov2c[(cov2c.pos<end_pos) & (cov2c.pos>start_pos)]
	# small_cov = small_cov.assign(pos=small_cov.pos-start_pos)
	seqs.append(small_seq)
	# covs.append(small_cov)

workers = multiprocessing.Pool(10)
tss_dist={}
for dist in workers.imap_unordered(kmer_distTSS, [ (seqs[i],args.k) for i in list(range(len(seqs)))]):
	for key,val in dist.items():
		if key not in tss_dist.keys():
			tss_dist[key]=[]
		for pos in val:
			tss_dist[key].append(pos)

oc_tbl=make_occurrencies_tbl(tss_dist)
# print "#hex"+'\t'.join([str(i) for i in list(oc_tbl)])
print("#hex"+'\t'.join([str(i) for i in list(oc_tbl)]))
for x in list(range(len(oc_tbl))):
	print(oc_tbl.index[x]+'\t'+'\t'.join([str(i) for i in oc_tbl.ix[x]]))
	# print oc_tbl.index[x]+'\t'+'\t'.join([str(i) for i in oc_tbl.ix[x]])	
















