import collections
import pysam as ps
import sys 
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count hexamers in BS converted fasta file.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('bam', type=str, help='Fasta input')
args = argparser.parse_args()

bamfile=args.bam
# bamfile='sorted_L3_trim1_R1_bismark_bt2_pe.deduplicated.bam'

trim_r1=9
trim_r2=8
bed=[]
with ps.AlignmentFile(bamfile,"rb") as bam:
	for r in bam.fetch(until_eof=True):
		if r.is_read1:
			a=r.qname
			bed.append((r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6))
		if r.is_read2:
			a=r.qname
			bed.append((r.reference_name, r.pos - trim_r2, r.pos - trim_r2 + 6))

for line in bed:
	print('\t'.join([str(i) for i in line]))
			