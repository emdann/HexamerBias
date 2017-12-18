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


def bestAlignmentKmers(params):
	seq,fqHex = params
	spl_seq=[['T','C'] if i == "C" else i for i in list(seq)]
	# probs=[[1-df.frac[df.pos==i+1].values[0], df.frac[df.pos==i+1].values[0]] if i in list(df["pos"]-1) else [0.99,0.01] if len(spl_seq[i])==2 else [1] for i in list(range(len(spl_seq)))]
	# bases=[[spl_seq[i],'C'] if i in list(df[(df.pos==i+1) & (df.strand=='+')].pos-1) else [spl_seq[i],'G'] if i in list(df[(df.pos==i+1) & (df.strand=='-')].pos-1) else spl_seq[i] for i in list(range(len(spl_seq)))]
	hex=spl_seq	
	hex_perm = list(it.product(*hex))
	perms=[''.join(i) for i in hex_perm]
	hex_perm.extend([tuple(str(Seq(seq, generic_dna).complement())) for seq in perms])
	maxScore=0
	for i in range(len(hex_perm)):
		al = pairwise2.align.globalxx(fqHex, ''.join(hex_perm[i]))
		score = al[0][2]
		if score > maxScore:
			maxScore = score
	return(fqHex, maxScore)

fasta='/hpc/hub_oudenaarden/edann/genomes/mm10/chr10.fa.gz'
bedfile='/hpc/hub_oudenaarden/edann/primeGen_L1.bed'

chromosome = fasta.split('/')[-1].split('.')[0]
bed = pd.read_csv(bedfile, sep="\t", header=None, dtype={4:int}, names=["chrom", "start", "end", "hex"])
refgen = bed[bed.chrom==chromosome]

# chrs = bed.chrom.unique()
# records = SeqIO.to_dict(SeqIO.parse(open(fasta), 'fasta'))

with ps.FastxFile(fasta) as chr:
 	for entry in chr:
 		seq=entry.sequence.upper()

hexs=list(refgen.hex)
covs=[]
seqs=[]
start = refgen.start
for i in start:
	start_pos = i
	end_pos = i+6
	small_seq=seq[start_pos:end_pos]
	seqs.append(small_seq)

workers = multiprocessing.Pool(10)
counts=[]

# track=1
for hex,bestScore in workers.imap_unordered(bestAlignmentKmers, [ (seqs[i],hexs[i]) for i in list(range(len(seqs)))]):
        # print("Adding count no. "+str(track))
        counts.append((hex,bestScore))
        # track+=1

for kmer, abundance in counts.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")
