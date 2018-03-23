import argparse
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Make pt counts tables for bs-seq \n By Emma Dann")
argparser.add_argument('abfile', type=str, help='Csv file of kmer abundance on reference genome')
argparser.add_argument('covfile', type=str, help='predicted coverage file')
argparser.add_argument('bed', type=str, help='bed of regions OI')
argparser.add_argument('refgen', type=str, help='Fasta file of reference genome')
args = argparser.parse_args()

abundanceFile = args.abfile
covFile = args.covfile
refgen = args.refgen
bedFile = args.bed

# Read files
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
coverage = pd.read_table(covFile, index_col=0, sep='\t',compression=findCompr(covFile))
with open(bedFile, 'r') as f:
    beds = [line.strip() for line in f.readlines()]

# Compute artificial coverage
density = template_density(coverage.exp,abundance)
covBed = artificial_cov(beds,refgen,density)
save_coverage_bed(covBed, outfile = bedFile.split('.bed')[0]+'.artCov.bed')
