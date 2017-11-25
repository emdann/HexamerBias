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
argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
argparser.add_argument('refgen', type=str, help='RefGen file of chr of interest (downloaded from UCSC genome browser)')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
args = argparser.parse_args()

def kmer_distTSS(params):
	seq,df,tss,k = params
	spl_seq=list(seq)
	bases=[[spl_seq[i],'T'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'A'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]	
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

chr=args.fasta.split('/')[-1].split('.')[0]
# print chr
refgen = pd.read_csv(args.refgen, sep="\t", usecols=[0,2,3,4], header=0, dtype={4:int}) 
refgen=refgen[refgen.chrom==chr]
refgen = refgen.drop_duplicates(subset=None, keep='first', inplace=False)
tss=refgen.txStart
# print refgen

with ps.FastxFile(args.fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

cov2c = pd.read_csv(args.cov2c, sep="\t", header=None) 
cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
cov2c=cov2c.assign(frac=cov2c.C / (cov2c.C+cov2c["T"]))
cov2c=cov2c[-np.isnan(cov2c.frac)]

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
# print "Covs computed,"+str(len(covs))+" tables"
# tss = refgen.txStart
# tss=np.array([random.randrange(1000) for x in range(100)])
# save_obj(tss_dist, TSS)

workers = multiprocessing.Pool(5)
tss_dist={}
for dist in workers.imap_unordered(kmer_distTSS, [ (seqs[i],covs[i],tss,args.k) for i in list(range(len(seqs)))]):
	for key,val in dist.items():
		if key not in tss_dist.keys():
			tss_dist[key]=[]
		for pos in val:
			tss_dist[key].append(pos)

save_obj(tss_dist, "TSS"+args.cov2c)

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

oc_tbl=make_occurrencies_tbl(tss_distance)
save_obj(oc_tbl, occurrencies_tbl)
















