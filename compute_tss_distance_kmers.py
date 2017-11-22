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

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('fasta', type=str, help='chr.fasta input')
argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
argparser.add_argument('refgen', type=str, help='RefGen file of chr of interest (downloaded from UCSC genome browser)')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

# refgen = pd.read_csv(args.refgen, sep="\t", usecols=[0,2,3,4],  nrows=100, header=0, dtype={4:int}) 

# refgen = pd.read_csv("~/mnt/edann/hexamers/hgTables.refgen.mm10.txt", sep="\t", usecols=[0,2,3,4],  nrows=100, header=0, dtype={4:int}) 
# refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)

# with ps.FastxFile("mnt/edann/genomes/mm10/chr1.fa") as chr:
#  	for entry in chr:
#  		seq=entry.sequence.upper()

# cov2c = pd.read_csv("mnt/edann/hexamers/test/test_chr1.txt", sep="\t", header=None) 
# cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
# cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
# df=cov2c[-np.isnan(cov2c.frac)]

# spl_seq=list(seq)
# bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]

# tss=np.array([random.randrange(1000) for x in range(100)])
# [[spl_seq[i],'T'] if i in list(df1000[df1000.pos==i+1]) for i in list(range(len(spl_seq1000)))]

def kmer_distTSS(bases,tss, k=6):
	tss_dist={}
	for i in range(0, len(bases)-k+1):
	    hex = bases[i:i+k]
	    hex_perm = list(it.product(*hex))
	    for x in range(len(hex_perm)):
	    	if(''.join(hex_perm[x]) not in tss_dist.keys()):
	    		tss_dist[''.join(hex_perm[x])]=[]
	    	tss_dist[''.join(hex_perm[x])].append(i-tss[np.argmin(abs(i-tss))])
	return(tss_dist)

# plt.hist(tss_dist["TTTTGT"])
# plt.show()

# # TO PLOT
# from scipy.stats import gaussian_kde
# density1 = gaussian_kde(tss_dist["TTTTGT"])
# density2 = gaussian_kde(tss_dist["TTTTTG"])
# xs = np.linspace(-10,10)
# density.covariance_factor = lambda : .00001
# density._compute_covariance()
# plt.plot(xs,density1(xs))
# plt.plot(xs,density2(xs))
# plt.show()

# # Make sparse matrix from dict 
# row_ind = [k for k, v in d.items() for _ in range(len(v))]
# col_ind = [i for ids in d.values() for i in ids]
# col_ind = list(set(col_ind))
# df=pd.DataFrame(np.zeros(shape=(len(row_ind),len(col_ind))))
# df.columns=col_ind
# df.index=row_ind
# for key,val in d.items():
# 	for i in val:
# 		print(key,i)

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def compute_tss_distance_kmers():
	refgen = pd.read_csv(args.refgen, sep="\t", usecols=[0,2,3,4],  nrows=100, header=0, dtype={4:int}) 
	refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)
	with ps.FastxFile(args.fasta) as chr:
	 	for entry in chr:
	 		seq=entry.sequence.upper()
	spl_seq=list(seq)
	cov2c = pd.read_csv(args.cov2c, sep="\t", header=None) 
	cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
	cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
	df=cov2c[-np.isnan(cov2c.frac)]
	bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]	
	# tss = refgen.txStart
	tss=np.array([random.randrange(1000) for x in range(100)])
	tss_dist = kmer_distTSS(bases,tss)
	# save_obj(tss_dist, TSS)
	return(tss_dist)

def make_occurrencies_tbl(tss_dist):
	dic_oc={}
	hexs=list(tss_dist.keys())
	counts=[collections.Counter(i) for i in tss_dist.values()]
	col_ind = [i for ids in tss_dist.values() for i in ids]
	col_ind = list(set(col_ind))
	for i in list(range(len(hexs))):
		dic_oc[hexs[i]]=counts[i]
	for key,val in dic_oc.items():
		for i in col_ind:
			if i not in list(val.keys()):
				val[i]=0
	oc_tbl=pd.DataFrame(dic_oc).T
	return(oc_tbl)

tss_distance=compute_tss_distance_kmers()
oc_tbl=make_occurrencies_tbl(tss_distance)
save_obj(oc_tbl, occurrencies_tbl)
















