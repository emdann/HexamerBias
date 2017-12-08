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

## rows= first nucleotide (AGTC) cols=second nucleotide (AGTC)

rand_hex=[''.join(i) for i in list(it.product(list("ATCG"),repeat=6))]

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

#rand_hex_deltaG={}
for hex in rand_hex:
	print(hex, compute_deltaG(hex)) 
	#rand_hex_deltaG[hex]=compute_deltaG(hex)
	



























