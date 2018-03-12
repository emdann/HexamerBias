import pysam as ps
import argparse
import os
import pybedtools as pbt

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Extract position of primer placement from trimmed section of aligned reads. By Emma Dann")
argparser.add_argument('bam', type=str, help='Input bam file')
argparser.add_argument('refgen', type=str, help='Fasta of reference genome (needs to be unzipped)')
argparser.add_argument('type', type=str, help='Original data of bam (bs_pe, bs_se or rna)')
argparser.add_argument('-o', type=str, help='path to directory to save output')
args = argparser.parse_args()

def get_template_bed(bamfile, type, trim=9):
	'''
	Extract positions of template regions for primers of aligned reads in bam.
	For type = 'rna' extracts the first 6 bases of the aligned read
	'''
	if type=='bs_pe':
		trim_r1=9
		trim_r2=8
		bed=[]
		with ps.AlignmentFile(bamfile,"rb") as bam:
			for r in bam.fetch(until_eof=True):
				if r.is_read1:
					bed.append((r.reference_name, r.pos - trim_r1, r.pos - trim_r1 + 6, r.qname, 1))
				# if r.is_read2:
		# 	bed.append((r.reference_name, r.pos - trim_r2, r.pos - trim_r2 + 6, r.qname, 2))
	if type=='bs_se':
		trim=10
		bed=[]
		with ps.AlignmentFile(bamfile,"rb") as bam:
			for r in bam.fetch(until_eof=True):
				bed.append((r.reference_name, r.pos - trim, r.pos - trim + 6, r.qname, 1))
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
	with open(outpath + outfile, 'w') as out:
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

bamfile = args.bam
type = args.type
fi = args.refgen
outpath = args.o

get_template_fasta(bamfile, fi, outpath, type)
