import pysam as ps
import argparse
import collections
import pandas as pd

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('fasta1', type=str, help='Fasta input')
argparser.add_argument('fasta2', type=str, help='Fasta input')
argparser.add_argument('bam', type=str, help='cytosine report input')
# argparser.add_argument('-t', type=int, default=8, required=False, help='Amount of threads to use.')
args = argparser.parse_args()

bam=args.bam
fasta1=args.fasta1
fasta2=args.fasta2

# bam='sorted_L3_trim1_R1_bismark_bt2_pe.deduplicated.bam'
# fasta1='L3_R1.fastq.gz'
# fasta2='L3_R2.fastq.gz'

ids_r1={}
ids_r2={}
with ps.AlignmentFile(bam,"rb") as bam:
	for r in bam.fetch(until_eof=True):
		if r.is_read1:
			a=r.qname
			ids_r1[a[:a.index('_')]]=r.flag
		else:
			a=r.qname
			ids_r2[a[:a.index('_')]]=r.flag

fq1_dict={99:[], 163:[], 147:[], 83:[]}
with ps.FastxFile(fasta1) as fastq:
 	for entry in fastq:
 		if entry.name in ids_r1.keys():
	 		fq1_dict[ids_r1[entry.name]].append(entry.sequence[:6])

fq2_dict={99:[], 163:[], 147:[], 83:[]}
with ps.FastxFile(fasta2) as fastq:
 	for entry in fastq:
 		if entry.name in ids_r2.keys():
	 		fq2_dict[ids_r2[entry.name]].append(entry.sequence[:6])

tot={99:[], 163:[], 147:[], 83:[]}
for dict in fq1_dict,fq2_dict:
	for k,v in dict.items():
		tot[k].extend(v)

# fq2_dict={}
# with ps.FastxFile(fasta2) as fastq:
#  	for entry in fastq:
#  		fq2_dict[entry.name]=entry.sequence[:6]


CTOT=pd.DataFrame.from_dict(collections.Counter(tot[147]), orient='index')
CTOB=pd.DataFrame.from_dict(collections.Counter(tot[163]), orient='index')
OT=pd.DataFrame.from_dict(collections.Counter(tot[99]), orient='index')
OB=pd.DataFrame.from_dict(collections.Counter(tot[83]), orient='index')
df=pd.concat([OT,OB,CTOB,CTOT], axis=1, join_axes=[OT.index])
df.columns=["OT", "OB", "CTOB", "CTOT"]

df.to_csv("/hpc/hub_oudenaarden/edann/"+fasta1.split('_')[0]+".splitHex")
# hex_count2=collections.Counter()
# for i in ids_r2:
# 	hex_count2[fq2_dict[i]] += 1 

# hex_count=hex_count1+hex_count2

# for kmer, abundance in hex_count1.most_common(): # sorts by abundance
# 	print(f"{kmer}\t{abundance}")
