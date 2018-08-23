import pysam as ps
import re
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
    description="Removes transcriptome contamination from raw reads from integrated Nla+CelSeq2 protocol. By Emma Dann")
argparser.add_argument('fastq', type=str, help='Fastq input file')
args = argparser.parse_args()

# fastqFile="CG-bulk100-hybES-NlaBsTx-NlaFraction-1_lmerged_R1.fastq.gz"
fastqFile=args.fastq
outname = fastqFile.split('.fastq')[0]+".noCelSeq.fastq.gz"

infile = ps.FastxFile(fastqFile, "r")
outfile = ps.FastxFile('ciao.fastq', 'w', template=infile)

for entry in infile:
    if "T"*15 not in entry.sequence:
        outfile.write(entry)
