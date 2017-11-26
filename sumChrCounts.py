import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import argparse
import scipy.sparse as sp
import fnmatch

dir='/home/emma/mnt/edann/hexamers'
# dir="/hpc/hub_oudenaarden/edann/hexamers"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, 'count_noNan_chr*'):
        files.append(dir+"/"+file)
# with open(dir+"/"+"files.txt",'rb') as f:
# 	for file in f.readlines():
# 		files.append(dir+"/"+file.strip("\n"))

hexAbundance=collections.Counter()
for file in files:
	with open(file,'rb') as f:
		for line in f:
			aline=line.strip('\n').split('\t')
			hexAbundance[aline[0]]+=float(aline[1])

for kmer, abundance in hexAbundance.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")
	# print kmer, abundance