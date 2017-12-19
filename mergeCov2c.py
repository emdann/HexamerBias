import pandas as pd
import sys 
import os
import argparse
import fnmatch
import collections
import gzip

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('chrom', type=str, help='chromosome cytosine report input')
args = argparser.parse_args()

chrom=args.chrom
dir='/hpc/hub_oudenaarden/edann/hexamers/kaester/met_extraction'
# dir="/hpc/hub_oudenaarden/edann/reikVSnla/met_extraction"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, 'cov2c*bismark.cov.gz'):
        files.append(dir+"/"+file)

cov_dict=collections.OrderedDict()
h=0
with open(files[0], "rb") as f:
	for line in f:
		# print("processing line", h)
		# h+=1
		if chrom+'\\t' in str(line):
			aline=line.decode().strip('\n').split('\t')
			cov_dict['\t'.join([aline[i] for i in [0,1,2,5,6]])]=[int(i) for i in aline[3:5]]
		
for file in files[1:]:
	with open(file,'rb') as f:
		for line in f:
			if chrom+'\\t' in str(line):
				aline=line.decode().strip('\n').split('\t')
				if '\t'.join([aline[i] for i in [0,1,2,5,6]]) not in cov_dict.keys():
					cov_dict['\t'.join([aline[i] for i in [0,1,2,5,6]])]=[0,0]
				cov_dict['\t'.join([aline[i] for i in [0,1,2,5,6]])][0]+= int(aline[3])
				cov_dict['\t'.join([aline[i] for i in [0,1,2,5,6]])][1]+= int(aline[4])

output_file=dir+'/merged_reference_CX.cov2c.'+chrom
with open(output_file, "w") as output:
	for key,val in cov_dict.items():
		list_key=key.split("\t")
		newLine =  list_key[:3]+[str(i) for i in val]+list_key[4:]
		print('\t'.join(newLine), file=output)

