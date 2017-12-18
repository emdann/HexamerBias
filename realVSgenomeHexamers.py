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

# bamfile=args.bam
# fasta1=args.fasta1
# fasta2=args.fasta2

bamfile='sorted_L3_trim1_R1_bismark_bt2_pe.deduplicated.bam'
fasta1='L3_R1.fastq.gz'
fasta2='L3_R2.fastq.gz'
ref='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'

records = SeqIO.to_dict(SeqIO.parse(open(ref), 'fasta'))
trim_r1=9
trim_r2=8
ids_r1={}
ids_r2={}
with ps.AlignmentFile(bamfile,"rb") as bam:
	for r in bam.fetch(until_eof=True):
		if r.is_read1:
			a=r.qname
			chr,start,end = (r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6)
			ids_r1[a[:a.index('_')]] = str(records[chr][start-1:end-1].seq).upper()
		if r.is_read2:
			a=r.qname
			chr,start,end = (r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6)
			ids_r2[a[:a.index('_')]] = str(records[chr][start-1:end-1].seq).upper()

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


