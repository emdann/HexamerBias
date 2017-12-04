import pandas as pd
import sys 
import os
import argparse
import fnmatch
import collections
import gzip

dir="/hpc/hub_oudenaarden/edann/hexamers"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, 'TSS_tab*'):
        files.append(dir+"/"+file)


tss_dist_dic={}
with open(files[0], "rb") as f:
	for line in f:
		aline=line.decode().strip('\n').split('\t')
		if len([float(i) for i in aline[1:]])==5995:
			tss_dist_dic[aline[0]]=[float(i) for i in aline[1:]]
		else:
			print("Error! Unexpected number of values in line "+ aline[0])
		# cov_dict['\t'.join([aline[i] for i in [0,1,2]])]=[int(i) for i in aline[4:]]

for file in files[1:]:
	with open(file,'rb') as f:
		for line in f:
			aline=line.decode().strip('\n').split('\t')
			if aline[0] not in tss_dist_dic.keys():
				tss_dist_dic[aline[0]]=[0]*5995
			if len([float(i) for i in aline[1:]])==5995:
				print("Processing hex "+aline[0]+"...")	
				tss_dist_dic[aline[0]]=list(pd.Series(tss_dist_dic[aline[0]])+pd.Series([float(i) for i in aline[1:]]))
			else:
				print("Error! Unexpected number of values in line "+ aline[0])

output_file=dir+"/sumTSS_distances_zf.txt"
with open(output_file, "w") as output:
	for key,val in tss_dist_dic.items():
		print(key+','+','.join([str(v) for v in val]), file=output)