import pysam as ps
import argparse
import collections

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('fasta1', type=str, help='Fasta input')
# argparser.add_argument('fasta2', type=str, help='Fasta input')
argparser.add_argument('bam', type=str, help='cytosine report input')
# argparser.add_argument('-t', type=int, default=8, required=False, help='Amount of threads to use.')
args = argparser.parse_args()

bam=args.bam
fasta=args.fasta1

# bam='/hpc/hub_oudenaarden/mauro/D50_D51_BS_2017_11/avo_files/D50-a-NN_tr2_R1_bismark_bt2.deduplicated.bam'
# fasta1='/hpc/hub_oudenaarden/mauro/D50_D51_BS_2017_11/
# fasta='./edann/crypts_bs/VAN1667/L1_R1.fastq.gz'

ids_r1=[]
# ids_r2=[]
with ps.AlignmentFile(bam,"rb") as bam:
	for r in bam.fetch(until_eof=True):
		# if r.is_read1:
		a=r.qname
		ids_r1.append(a[:a.index('_')])
		# else:
		# 	a=r.qname
		# 	ids_r2.append(a[:a.index('_')])

fq1_dict={}
with ps.FastxFile(fasta1) as fastq:
 	for entry in fastq:
 		fq1_dict[entry.name]=entry.sequence[:6]

# fq2_dict={}
# with ps.FastxFile(args.fasta2) as fastq:
#  	for entry in fastq:
#  		fq2_dict[entry.name]=entry.sequence[:6]

hex_count1=collections.Counter()
for i in ids_r1:
	hex_count1[fq1_dict[i]] += 1 

# hex_count2=collections.Counter()
# for i in ids_r2:
# 	hex_count2[fq2_dict[i]] += 1 

# hex_count=hex_count1+hex_count2

for kmer, abundance in hex_count1.most_common(): # sorts by abundance
	print(f"{kmer}\t{abundance}")
