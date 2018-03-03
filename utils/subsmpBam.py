import pysam as ps
import random
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Subsample bam file\n By Emma Dann")
argparser.add_argument('bamfile', type=str, help='bam input')
argparser.add_argument('bamout', help='bamout')
argparser.add_argument('-f', type=float, default=0.001, required=False, help='Subsample fraction size')
# argparser.add_argument('--chr', type=str, required=False, help='chromosome')
args = argparser.parse_args()

bamfile=args.bamfile
bamout=args.bamout
fraction=args.f

bam = ps.AlignmentFile(bamfile, 'rb') # Change me
output = ps.AlignmentFile(bamout, "wb", template=bam) # Change me

for read in bam.fetch(until_eof=True):
    if random.random() < fraction:
        output.write(read)

bam.close()
output.close()
