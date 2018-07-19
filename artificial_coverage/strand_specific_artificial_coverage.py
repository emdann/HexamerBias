import argparse
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

argparser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Build predicted genomic coverage track in regions of interest given a certain binding fraction for each sequence \n By Emma Dann",
    epilog='''\
        THE BED REGIONS MUSTN'T OVERLAP! Otherwise the entries won't be added in the right order and the
        programme will crash right at the end./)

        You can make sure there is no overlap running:
        bedtools sort -i myregions.bed | bedtools spacing -i stdin | awk '$4!=0'

    ''')
argparser.add_argument('abfile', type=str, help='Csv file of kmer abundance on reference genome')
argparser.add_argument('covfile', type=str, help='predicted coverage file (in .csv, template sequences in first column)')
argparser.add_argument('bed', type=str, help='bed of regions OI (NON OVERLAPPING INTERVALS!!!)')
argparser.add_argument('refgen', type=str, help='Fasta file of reference genome')
argparser.add_argument('--BS', type=str, default='no',help='BS conversion mode')
argparser.add_argument('--output', type=str, default='bigWig',help='Format of output file: bedGraph or bigWig')
argparser.add_argument('-t', type=int, default=10, required=False, help='Amount of threads to use.')
# argparser.add_argument('--numReads', type=int, default=10, required=False, help='Amount of threads to use.')
args = argparser.parse_args()

def save_bw_read_extend(beds,refgen_fasta,density, outfile,bs='no',readLength=10,threads=10):
    workers = multiprocessing.Pool(threads)
    bw = pbw.open(outfile, 'w')
    header=make_BigWig_header(refgen_fasta)
    chroms = [e[0] for e in header]
    bw.addHeader(header)
    intervals = []
    for chrom,posDic in workers.imap_unordered(artificial_cov_bed_entry, [(bed, refgen_fasta, density, readLength,bs) for bed in beds]):
        intervals.extend([(chrom,pos,val) for pos,val in kernel_smoothing(posDic, sigma=20).items()])
    print('--- Saving BigWig ---', flush=True)
    for chrom in chroms:
        print('Writing entries from '+ chrom, flush=True)
        exChr = sorted([i for i in intervals if i[0]==chrom])
        if len(exChr)>0:
            bw.addEntries(chrom, [k[1] for k in exChr], values=[k[2] for k in exChr], span=1)
    bw.close()
    return(intervals)

abundanceFile = args.abfile
covFile = args.covfile
refgen = args.refgen
bedFile = args.bed
outtype = args.output
bs=args.BS

# Read files
print('--- Reading input files ---', flush=True)
abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
coverage = pd.read_csv(covFile, index_col='template',compression=findCompr(covFile))
with open(bedFile, 'r') as f:
    beds = [line.strip() for line in f.readlines()]

# Compute artificial coverage
print('--- Computing density ---', flush=True)
density = template_density(coverage.exp,abundance)
save_bw_read_extend(beds,refgen,density, bedFile.split('.bed')[0] + covFile.split('.csv')[0] + '.artCov.bw',bs=bs, readLength=75, threads=args.t)
