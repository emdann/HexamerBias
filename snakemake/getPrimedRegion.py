import pysam as ps
import argparse
import os
import pybedtools as pbt

def get_template_bed(bamfile, type):
	'''
	Extract positions of template regions for primers of aligned reads in bam.
	For type = 'rna' extracts the first 6 bases of the aligned read
	'''
	if type=='bs':
		trim_r1=9
		trim_r2=8
		bed=[]
		with ps.AlignmentFile(bamfile,"rb") as bam:
			for r in bam.fetch(until_eof=True):
				if r.is_read1:
					bed.append((r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6, r.qname, 1))
				# if r.is_read2:
				# 	bed.append((r.reference_name, r.pos - trim_r2, r.pos - trim_r2 + 6, r.qname, 2))
	elif type=='rna':
		bed=[]
		with ps.AlignmentFile(bamfile,"rb") as bam:
			for r in bam.fetch(until_eof=True):
				if r.flag==0:
					bed.append((r.reference_name, r.pos, r.pos + 6, r.qname))
	return(bed)

def save_bedfile(bedlist, bamfile, outpath):
	'''
	Saves bedfile to defined output path
	'''
	sample = bamfile.split('/')[-1].split('.')[0]
	outfile = sample + '.primedreg.bed'
	with open(outpath +'/'+ outfile, 'w') as out:
		for line in bedlist:
			print('\t'.join([str(i) for i in line]), file=out)
	return(outpath + outfile)

def get_template_fasta(bamfile, fi, outpath, type):
	'''
	Makes fasta file of template sequences from bed file.
	Saves output fasta file.
	'''
	sample = bamfile.split('/')[-1].split('.')[0]
	bedfile = outpath + sample + '.primedreg.bed'
	if not os.path.exists(bedfile):
		bed = get_template_bed(bamfile, type=type)
		bedfile = save_bedfile(bed, bamfile, outpath)
	bed = pbt.BedTool(bedfile)
	faout = bed.sequence(fi, name=True)
	faout.save_seqs(bedfile.split('.bed')[0] + '.fa')
	return(faout)
