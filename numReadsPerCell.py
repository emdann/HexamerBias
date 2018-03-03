import pandas as pd
import pysam as ps
import collections
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Counts number of reads per cell in scRNA-seq alignment file (from barcode) \n By Emma Dann")
argparser.add_argument('bam', type=str, help='Bam input')
argparser.add_argument('-o', type=str, required=False, help='Output directory')
args = argparser.parse_args()

bamfile = args.bam
outpath = args.o

def num_reads_per_cell(bamfile, to_dir = None):
    '''
    Makes dictionary of number of aligned reads with high quality hexamer
    sequencing for each cell
    '''
    reads = []
    with ps.AlignmentFile(bamfile,"rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.flag==0 and not any(q<32 for q in r.query_qualities[:6]):
                reads.append(r.qname)
    cells = [r.split(':')[-1] for r in reads]
    numReads = collections.Counter(cells)
    if to_dir:
        sample = bamfile.split('/')[-1].split('.')[0]
        outfile = to_dir + '/' + sample + '.numReads.txt'
        with open(outfile, 'w') as f:
            print('cell', 'numReads', sep='\t', file = f)
            for cell,num in numReads.items():
                print(cell, num, sep='\t', file = f)
    return(numReads)

num_reads_per_cell(bamfile, to_dir = outpath)
