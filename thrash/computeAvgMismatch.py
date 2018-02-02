from hexVSprimed import *
import collections
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute average mismatch per hexamer type\n By Emma Dann")
argparser.add_argument('mmfile', type=str, help='Fasta input')
argparser.add_argument('--bs', help='Tolerate BS specific mismatches (C->T, G->A)')
args = argparser.parse_args()

mmfile=args.mmfile

if args.bs:
    tolerateBS=True
    outputfile=mmfile.replace('.mismatch.txt','')+'_avgMM_tolerateBSmm.txt'
else:
    tolerateBS=False
    outputfile=mmfile.replace('.mismatch.txt','')+'_avgMM.txt'

newDic={}
with open(mmfile, 'rt') as f:
    f.readline()
    for line in f.readlines():
        newDic[line.split()[0]]=[line.split()[1], line.split()[2]]

kmerMM={}
for primedSeq,hex in newDic.values():
    if hex not in kmerMM.keys():
        kmerMM[hex]=[]
    kmerMM[hex].append(mismatchKmers((primedSeq,hex), tolerateBSmm=tolerateBS))

with open(outputfile,'w') as f:
    for k,i in kmerMM.items():
        print(k,np.mean(i), sep='\t', file=f)
