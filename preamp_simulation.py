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
import random
from freeEnergy import compute_deltaG

seq='ATCGTTTACCG'
cov2c=pd.DataFrame({"pos":[3,10], "frac":[0.2,0.9], "strand":["+","+"]})
hexAb=strandSpecificCount((seq, cov2c,6))


def compute_prob_capture(h, temp=277):
	deltaG=compute_deltaG(h)
	N=hexAb[h]
	frac=N/sum(hexAb.values())
	pCapture=np.exp((1/(0.001986*(273+temp))*deltaG))
	return(pCapture)

# randProb={}
# hexs=[seq[i:i+6] for i in range(len(seq)-6+1)]
# for hex in set(hexs):
# 	randProb[hex]=random.random()

def amplification(hexs,it=1000):
	amp1=[]
	for h in hexs:
		# p=compute_prob_capture(h)
		p=np.exp(-compute_deltaG(h))*0.0001
		t=random.random()
		if t<p:
			amp1.append(h)
	return(amp1)

revcSeq=str(Seq(seq, generic_dna).reverse_complement())

def preamp_simulation(seq,revcSeq, rounds=5, it=1000):
	hexs=[seq[i:i+6] for i in range(len(seq)-6+1)]
	revHexs=[revcSeq[i:i+6] for i in range(len(revcSeq)-6+1)]
	k=0
	while k<it:
		amp1=amplification(hexs,it)+amplification(revHexs,it)
		r=2
		while r < rounds: 
			amp=amplification(amp1)
			amp1=amp
			r+=1
		k+=1
	ampMean=collections.Counter()
	for key,val in collections.Counter(amp1).items():
		ampMean[key]=float(val)/float(k)
	return(ampMean)



