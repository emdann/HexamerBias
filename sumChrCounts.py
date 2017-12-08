import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import scipy.sparse as sp
import fnmatch

dir="/hpc/hub_oudenaarden/edann/hexamers"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, 'count_convRef.chr*'):
        files.append(dir+"/"+file)

hexAbundance=collections.Counter()
for file in files:
with open(file,'rb') as f:
	for line in f:
		aline=line.decode().strip('\n').split('\t')
		hexAbundance[aline[0]]+=float(aline[1])

output_file=dir+"/rand_hex_abundance_with5amp.txt"
with open(output_file, "w") as output:
	for kmer, abundance in hexAbundance.most_common(): # sorts by abundance
		print(f"{kmer}\t{abundance}", file=output)
	# print kmer, abundance