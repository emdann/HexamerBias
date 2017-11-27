import pandas as pd
import sys 
import os
import argparse
import fnmatch
import collections
import gzip

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
#argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
args = argparser.parse_args()

dir='/home/emma/mnt/edann/reikVSnla/met_extraction'
# dir="/hpc/hub_oudenaarden/edann/reikVSnla/met_extraction"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, '*.cov.gz'):
        files.append(dir+"/"+file)

cov_dict=collections.OrderedDict()
with gzip.open(files[0], "rb") as f:
	for line in f:
		aline=line.decode().strip('\n').split('\t')
		cov_dict['\t'.join([aline[i] for i in [0,1,2]])]=[int(i) for i in aline[4:5]]

for file in files[1:]:
	with gzip.open(file,'rb') as f:
		for line in f:
			aline=line.decode().strip('\n').split('\t')
			cov_dict['\t'.join([aline[i] for i in [0,1,2]])][0]+= int(aline[4])
			cov_dict['\t'.join([aline[i] for i in [0,1,2]])][1]+= int(aline[5])

for key,val in cov_dict.items():
	list_key=key.split("\t")
	newLine =  list_key[:3]+[str(i/2) for i in val]+list_key[3:]
	print('\t'.join(newLine))


