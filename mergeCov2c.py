import pandas as pd
import sys 
import os
import argparse
import fnmatch

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute distance to TSS of hexamers in BS converted chromosome.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('cov2c', type=str, help='chromosome cytosine report input')
args = argparser.parse_args()

dir="/hpc/hub_oudenaarden/edann/hexamers/kaester/met_extraction"
files=[]
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, '*.CpG_report.txt.gz'):
        files.append(file)


cov2c = pd.read_csv(dir+"/"+files[0], sep="\t", header=None) 
cov2c.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]

for file in files[1:]:
	cov2c_2 = pd.read_csv(dir+"/"+file, sep="\t", header=None) 
	cov2c_2.columns = ["chr", "pos", "strand", "C", "T", "context", "flank"]
	## gives MemoryError --> ask Buys
	df=pd.merge(cov2c,cov2c_2,how="outer",on=["chr","pos","strand","context", "flank"] ).fillna(0)
	df=df.assign(C=df.C_x+df.C_y, T=df.T_x+df.T_y)
	cov2c=df[["chr","pos","strand","C","T","context","flank"]]