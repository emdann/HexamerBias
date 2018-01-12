import pandas as pd
import pysam as ps
import collections
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
args = argparser.parse_args()

def mismatchKmers(params):
    seq,fqHex = params
    a = pairwise2.align.globalxs(fqHex, seq, -1000, -1000) # alignment funciton to not open gaps
    score = a[0][2]
    mmatch=6-score
    return(mmatch)

def make_occurrencies_tbl(seqDict):
    dic_oc={}
    for seqs in seqDict.values():
        primed,hex = seqs[0],seqs[1]
        if hex not in dic_oc.keys():
            dic_oc[hex]=[]
        dic_oc[hex].append(primed)
    hexs=list(dic_oc.keys())
    countOcc={}
    counts=[collections.Counter(i) for i in dic_oc.values()]
    for i in list(range(len(hexs))):
        countOcc[hexs[i]]=counts[i]
    oc_tbl=pd.DataFrame(countOcc)
    oc_tbl=oc_tbl.fillna(0)
    return(oc_tbl)

fastqFile=args.fasta
primedregFile=args.primedreg

#fastqFile='/hpc/hub_oudenaarden/edann/crypts_bs/VAN1667/L1_R1.fastq.gz'
#primedregFile='/hpc/hub_oudenaarden/edann/hexamers/L1_R1_primed_seq.original.fa'

dict={}
with ps.FastxFile(primedregFile) as fastq:
    for entry in fastq:
        dict[entry.name.split('_')[0]]=[entry.sequence.upper()]

with ps.FastxFile(fastqFile) as fastq:
    for entry in fastq:
#        print(entry.name)
        if entry.name in dict.keys():
            dict[entry.name].append(entry.sequence[:6])

with open(primedregFile+'.mismatch.txt', 'w') as output:
    print('read', 'primedSeq', 'hex', sep='\t', file=output)
    for key,val in newdict.items():
        print(key, val[0], val[1], sep='\t', file=output)


tbl = make_occurrencies_tbl(dict)
tbl.to_csv(primedregFile+'.csv')
