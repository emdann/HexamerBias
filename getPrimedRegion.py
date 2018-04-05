print('whatever2')

import pysam as ps
import argparse
import os
import pybedtools as pbt

print('whatever')

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Extract position of primer placement from trimmed section of aligned reads. By Emma Dann")
argparser.add_argument('bam', type=str, help='Input bam file')
argparser.add_argument('refgen', type=str, help='Fasta of reference genome (needs to be unzipped)')
argparser.add_argument('-t', type=str, default='bs_pe',required=False, help='Original data of bam (bs_pe, bs_se or rna)')
argparser.add_argument('-s', action='store_true', help='require strandedness')
argparser.add_argument('-o', type=str, help='path to directory to save output')
args = argparser.parse_args()

print('whatever1')

## For some reason doesn't complete the fasta file when running everything in one.
# Consider splitting into two scripts and retesting

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
				bed.append((r.reference_name, r.pos+1 - trim, r.pos+1 - trim + 6, r.qname, 1))
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
	Makes fasta file of template sequences from bed file. Calling bedtools from terminal.
	Saves output fasta file.
	'''
	sample = bamfile.split('/')[-1].split('.')[0]
	bedfile = outpath + sample + '.primedreg.bed'
	if not os.path.exists(bedfile):
		print("Saving bed...")
		bed = get_template_bed(bamfile, type=type)
		bedfile = save_bedfile(bed, bamfile, outpath)
		print("Bed saved!")
	command = "/hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -fi " + fi + " -bed " + bedfile + " -fo " + bedfile.split('.bed')[0] + '.fa'
	print("Saving fasta...")
	os.system(command)
	print("fasta saved!")
	return('')


def get_template_fasta_pybedtools(bamfile, fi, outpath, type):
	'''
	Makes fasta file of template sequences from bed file.
	Saves output fasta file.
	'''
	## For some reason doesn't complete the fasta file when running everything in one.
	# Consider splitting into two scripts and retesting
	sample = bamfile.split('/')[-1].split('.')[0]
	bedfile = outpath + sample + '.primedreg.bed'
	if not os.path.exists(bedfile):
		bed = get_template_bed(bamfile, type=type)
		bedfile = save_bedfile(bed, bamfile, outpath)
	bed = pbt.BedTool(bedfile)
	faout = bed.sequence(fi, name=True)
	faout.save_seqs(bedfile.split('.bed')[0] + '.fa')
	return(faout)

## ---- STRAND SPECIFIC TEMPLATES ----
def get_strandspecific_template_bed(bamfile, trim=9, type='bs_se'):
	'''
	Extract positions of template regions for primers of aligned reads in bam (se mapped bismark).
	For type = 'rna' extracts the first 6 bases of the aligned read
	'''
	trim=trim
	bed=[]
	if type=='bs_se': # Define flags
		plus=0
		minus=16
	if type=='bs_pe':
		plus=99
		minus=83
	with ps.AlignmentFile(bamfile,"rb") as bam:
		for r in bam.fetch(until_eof=True):
			if r.flag==plus:
				bed.append((r.reference_name, r.pos+1 - trim, r.pos+1 - trim + 6, r.qname, '.', '+'))
			if r.flag==minus:
				bed.append((r.reference_name, r.pos+1 - trim, r.pos+1 - trim + 6, r.qname, '.', '-'))
	return(bed)

def get_strandspecific_template_fasta(bamfile, fi, outpath, type):
	'''
	Makes fasta file of template sequences from bed file. Calling bedtools from terminal.
	Saves output fasta file.
	'''
	sample = bamfile.split('/')[-1].split('.')[0]
	bedfile = outpath + sample + '.primedreg.bed'
	if not os.path.exists(bedfile):
		print("Saving bed...")
		bed = get_strandspecific_template_bed(bamfile, type=type)
		bedfile = save_bedfile(bed, bamfile, outpath)
		print("Bed saved!")
	command = "/hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -s -fi " + fi + " -bed " + bedfile + " -fo " + bedfile.split('.bed')[0] + '.fa'
	print("Saving fasta...")
	os.system(command)
	print("fasta saved!")
	return('')

print('again')

bamfile = args.bam
type = args.t
fi = args.refgen
outpath = args.o
strandedness = args.s
print(strandedness)

if strandedness:
	print('Finding strand specific templates...')
	get_strandspecific_template_fasta(bamfile, fi, outpath, type)
else:
	print('Finding non strand specific templates...')
	get_template_fasta(bamfile, fi, outpath, type)
