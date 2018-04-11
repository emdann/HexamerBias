import argparse
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts tables for bs-seq \n By Emma Dann")
argparser.add_argument('abfile', type=str, help='Csv file of kmer abundance on reference genome')
argparser.add_argument('covfile', type=str, help='predicted coverage file (in .csv, template sequences in first column)')
argparser.add_argument('bed', type=str, help='bed of regions OI')
argparser.add_argument('refgen', type=str, help='Fasta file of reference genome')
argparser.add_argument('--output', type=str, default='bigWig',help='Format of output file: bedGraph or bigWig')
argparser.add_argument('-t', type=int, default=10, required=False, help='Amount of threads to use.')
args = argparser.parse_args()

def artificial_cov_bed_entry(params):
    bedEntry,refgen_fasta,density,read_length=params
    chr,start,end = bedEntry.split()
    print('Processing entry ', bedEntry, flush=True)
    seq = ps.FastaFile(refgen_fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
    strandSpecificPosDic = sum_strands_per_base_cov(seq, density, int(start), readLength=read_length)
    smoothPosDic = kernel_smoothing(strandSpecificPosDic)
    return(chr,smoothPosDic)

def save_bw_read_extend(beds,refgen_fasta,density, outfile,readLength=10,threads=10):
    workers = multiprocessing.Pool(threads)
    bw = pbw.open(outfile, 'w')
    bw.addHeader(make_BigWig_header(refgen_fasta))
    for chrom,posDic in workers.map(artificial_cov_bed_entry, [(bed, refgen_fasta, density, readLength) for bed in beds]):
        # print(chrom)
        bw.addEntries(chrom, [k for k in posDic.keys()], values=[v for v in posDic.values()], span=1)
    bw.close()

abundanceFile = args.abfile
covFile = args.covfile
refgen = args.refgen
bedFile = args.bed
outtype = args.output

# Read files
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
coverage = pd.read_csv(covFile, index_col='template',compression=findCompr(covFile))
with open(bedFile, 'r') as f:
    beds = [line.strip() for line in f.readlines()]

# Compute artificial coverage
density = template_density(coverage.exp,abundance)
save_bw_read_extend(beds,refgen,density,bedFile.split('.bed')[0]+'.artCov.bw', readLength=75, threads=args.t)
