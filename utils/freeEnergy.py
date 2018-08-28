import itertools as it
import pandas as pd
import collections
import numpy as np
import sys 
import os
import math

## rows= first nucleotide (AGTC) cols=second nucleotide (AGTC)

# rand_hex=[''.join(i) for i in list(it.product(list("ATCG"),repeat=6))]

def compute_deltaG(hex):
	deltaG_tbl=np.array([
		-1.28,-1.54,-1.12,-1.72,
		-1.58,-2.07,-1.72,-2.53, 
		-0.85, -1.73,-1.28,-1.58, 
		-1.74,-2.49,-1.54,-2.07
		]).reshape((4,4))
	pos={"A":0,"G":1,"T":2,"C":3}
	hex_G=0
	for i in list(range(len(hex)-1)):
		hex_G += deltaG_tbl[pos[hex[i]], pos[hex[i+1]]]
	return(hex_G)

def compute_deltaG_ions(hex, NaConc=0.05, MgConc=0.01):
	deltaGNa_tbl=np.array([
		-1.24,-1.33,-1.07,-1.49,
		-1.52,-1.92,-1.49,-2.38, 
		-0.84, -1.76,-1.24,-1.52, 
		-1.76,-2.34,-1.33,-1.92
		]).reshape((4,4))
	mNa_tbl=np.array([
		0.142 , 0.087 , 0.107 , 0.115 ,
		0.142 , 0.069 , 0.115 , 0.080 , 
		0.089 , 0.096 , 0.142 , 0.142 , 
		0.096 , 0.117 , 0.087 , 0.069
		]).reshape((4,4))
	mMg_tbl=np.array([
		0.086 , 0.070 , 0.092 , 0.073 ,
		0.075 , 0.032 , 0.073 , 0.060 , 
		0.087 , 0.078 , 0.086 , 0.075 , 
		0.078 , 0.057 , 0.070 , 0.032
		]).reshape((4,4))
	deltaGMg_tbl=np.array([
		-1.69,-1.81,-1.55,-1.91,
		-1.88,-2.18,-1.91,-2.74, 
		-1.38, -2.17,-1.69,-1.88, 
		-2.17,-2.65,-1.81,-2.18
		]).reshape((4,4))
	mean_tbl = np.add(deltaGMg_tbl, deltaGNa_tbl)/2
	
	pos={"A":0,"G":1,"T":2,"C":3}
	hex_G=0
	for i in list(range(len(hex)-1)):
		hex_G += mean_tbl[pos[hex[i]], pos[hex[i+1]]] - mNa_tbl[pos[hex[i]], pos[hex[i+1]]]*math.log(NaConc) - mMg_tbl[pos[hex[i]], pos[hex[i+1]]]*math.log(MgConc)
	if hex[0] in ["A","T"]:
		hex_G = hex_G + 0.43
	if hex[-1] in ["A","T"]:
		hex_G = hex_G + 0.43
	return(hex_G)


# #rand_hex_deltaG={}
# for hex in rand_hex:
# 	print(hex, compute_deltaG_ions(hex)) 
# 	#rand_hex_deltaG[hex]=compute_deltaG(hex)
	



























