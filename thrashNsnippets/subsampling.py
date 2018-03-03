import pandas as pd
import random
import gzip
import collections
import numpy as np


file="/hpc/hub_oudenaarden/edann/reikVSnla/met_extraction/CpG_context_K562_nla_cat_R2_bismark_bt2.deduplicated.txt.gz"
output_file="/hpc/hub_oudenaarden/edann/reikVSnla/Nla_lorentz.txt"
# CpG_context = pd.read_csv(file,  compression='gzip', sep="\t", skiprows=1, usecols=[0], names=["reads"])

# reads=set(list(CpG_context.reads))


bismark_dic={}
with gzip.open(file, "rb") as f:
	for line in f:
		aline=line.decode().strip('\n').split('\t')
		if aline[0] not in bismark_dic.keys():
			bismark_dic[aline[0]]=[]
		bismark_dic[aline[0]].append('.'.join(aline[2:4]))
	f.close()

reads=list(bismark_dic.keys())
n=len(bismark_dic.keys())
subsmp_dic={}
for p in np.arange(1,4000000, 50000):
	counts=collections.Counter()
	sample = random.sample(reads, p)
	for i in sample:
		for k in bismark_dic[i]:
			counts[k]+=1
	subsmp_dic[str(p)]=len(counts.keys())

with open(output_file,"w") as output:
	for key,val in subsmp_dic.items():
		print(key+"\t"+str(val), file=output)
