from hexVSprimed import *
import pysam as ps
import collections
import argparse
import pandas as pd

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
args = argparser.parse_args()

fasta = args.primedseq_fasta

bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/zebrafish/gk2a-2.sam'
fasta='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2_primed_seq.fa'

templDic={}
with ps.FastxFile(fasta) as templ:
    for entry in templ:
        seq, name = entry.sequence.upper(), entry.name
        templDic[name]=[seq]


with ps.AlignmentFile(bamfile,"rb") as bam:
    for r in bam.fetch(until_eof=True): 
        if r.flag==0:
            templDic[r.qname].append(r.seq[0:6])
        elif r.flag==16:
            templDic.pop(r.qname)

SOI=["GTGTGT", "TGTGTG" ,"TTCCCC", "TAGAAC" ,"ACGAAT" ,"TAGGAA" ,
    "AAAAAA" ,"TTTTTT", "GGGGGG", "CCCCCC",
    "GTAGGT","GTACGC", "GTACCG", "AAAAAT", "ATGAAA",
    "CGACGT", "ATTCAC", "TGATTT" ,"CTTTAT" ,"ACGACT", "AACCCT", "GGATGC", "CTGGAT", "CAGCAG", "TCGCAC"]

SOICounter = {key:collections.Counter() for key in SOI}
for seqs in templDic.values():
    templ,primer=seqs[0],seqs[1]
    if primer in SOI:
        SOICounter[primer][templ]+=1

pwms={k:makePWMtempl(v) for k,v in SOICounter.items()}
entropies = [computeEntropy(i) for i in pwms]

for k,pwm in pwms.items():
    pwm.to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/'+k+'.primerDanRer.pwm.csv')

tbl = pd.DataFrame(SOICounter).fillna(0)
tbl.to_csv(primedregFile.split("/")[-1].split('.fa')[0] +'.csv')
