import json
import os
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Makes json config file for mapping snakemake. Has to be run from ")
# argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('-o', type=str, default='.', required=False, help='directory of fastq files (default: cwd)')
argparser.add_argument('--clip5prime', type=str, default='9', required=False, help='number of bases to clip from 5 prime end (default: 9)')
args = argparser.parse_args()

inputDir=args.o
clip=args.clip5prime

samples = list(set([f.split('_R')[0] for f in os.listdir(inputDir) if "fastq.gz" in f]))

configDic = {}
configDic["sample"]=samples
configDic["clip5prime"] = clip

with open('configfile.json', 'w') as outfile:
    json.dump(configDic, outfile)
