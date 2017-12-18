import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys 
import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count hexamers in BS converted fasta file.\n Do it per chromosome! By Emma Dann")
argparser.add_argument('bam', type=str, help='Bam input')
argparser.add_argument('fasta1', type=str, help='Bam input')
argparser.add_argument('fasta2', type=str, help='Bam input')
args = argparser.parse_args()

bamfile=args.bam
fasta1=args.fasta1
fasta2=args.fasta2

trim_r1=9
trim_r2=8
ids_r1={}
ids_r2={}
with ps.AlignmentFile(bamfile,"rb") as bam:
	for r in bam.fetch(until_eof=True):
		if r.is_read1:
			a=r.qname
			ids_r1[a[:a.index('_')]] = (r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6)
		if r.is_read2:
			a=r.qname
			ids_r2[a[:a.index('_')]] = (r.reference_name, r.pos - trim_r2, r.pos - trim_r2 + 6)

def fastqVSgenHexamers(fastqfile, ids):
	fq_dict = {}
	with ps.FastxFile(fastqfile) as fastq:
		for rec in fastq:
			if rec.name in ids.keys():
				fq_dict[rec.name]=rec.sequence[:6]
	for id,hex in ids.items():
		ids[id]=(hex, fq_dict[id])
	return(ids)

r1 = fastqVSgenHexamers(fasta1, ids_r1)
r2 = fastqVSgenHexamers(fasta2, ids_r2)

for line in r1.values():
	new_line=[i for i in line[0]]
	new_line.append(line[1])
	print("\t".join([str(i) for i in new_line]))
