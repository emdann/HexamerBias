import pandas as pd
import pysam as ps
import collections
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import math
import random

def mismatchKmers(params, tolerateBSmm):
    '''
    calculates the mismatch between the reference aligned sequence (seq) and the hexamer found in the fastq file (fqHex)
    tolerateBSoption to ignore C to T and G to A mismatch.
    '''
    seq,fqHex = params
    if tolerateBSmm==True:
        score=0
        for i in range(len(fqHex)):
            x,y=fqHex[i],seq[i]
            print(x,y)
            if x==y:
                score+=1
            elif x == 'T' and y =='C':
                score+=1
            elif x == 'G' and y =='A':
                score+=1
            else:
                score+=0
    else:
        a = pairwise2.align.globalxs(fqHex, seq, -1000, -1000) # alignment funciton to not open gaps
        score = a[0][2]
    mmatch=6-score
    return(mmatch)

def mmPerPosition(openbam):
    a = []
    for r in openbam.fetch(until_eof=True):
        if r.flag!=4:
            for i in r.get_aligned_pairs(with_seq=True):
                rpos,refbase=i[0],i[2]
                if refbase and refbase.islower():
                    mm=refbase.upper()+'->'+r.seq[rpos]
                    if r.flag==16:
                        strand='-'
                    if r.flag==0:
                        strand='+'
                    a.append({'type':mm, 'pos':rpos, 'strand':strand})
    return(a)

def mmPerQual(openbam):
    mmQual = {}
    for r in openbam.fetch(until_eof=True):
        if r.flag==0:
            for i in r.get_aligned_pairs(with_seq=True):
                rpos,refbase=i[0],i[2]
                if rpos and refbase:
                    mm=refbase.upper()+'->'+r.seq[rpos]
                    qual=r.query_qualities[rpos]
                    if mm not in mmQual.keys():
                        mmQual[mm]=[]
                    mmQual[mm].append(qual)
    return(mmQual)

def usageMmPerHex(openbam, ref, templDic=None):
    '''
    Computes usage and average number of mismatches and count of mismatching events for a specific hexamer
    either as primer or as template sequence
    (from RNAseq bwa aligned bam)
    '''
    mmHex={}
    for r in openbam.fetch(until_eof=True):
        if r.flag==0:
            refbases = [r.get_aligned_pairs(with_seq=True)[i][2] if r.query_qualities[i] >=32 else 'lowqual' for i in range(0,len(r.get_aligned_pairs())) ]
            seq = r.seq
            start,end = 0,6
            if ref=='primer':
                hex=seq[0:6]
                print('Primer seq: '+hex)
            elif ref=='template':
                hex = templDic[r.qname]
                print('Template seq: '+hex)
            mm = 0
            for b in refbases[start:end]:
                if not b:
                    mm+=1
                elif b!='lowqual' and b.islower():
                    mm+=1
            if hex not in mmHex.keys():
                mmHex[hex]=[]
            mmHex[hex].append(mm)
    hexDic={}
    for k,i in mmHex.items():
        mmEvents=[a for a in i if a!=0]
        hexDic[k]=(len(i),len(mmEvents), np.mean(i))
    return(hexDic)

def countHex(openbam):
    hexCounts=collections.Counter()
    for r in openbam.fetch(until_eof=True):
        if r.flag==0:
                seq = r.seq
                hex=str(seq)[0:6]
                hexCounts[hex]+=1
    return(hexCounts)

def make_mm_table(listDic):
    pos = [int(i['pos']) for i in listDic]
    occ = dict.fromkeys(list(range(0,max(pos)+1)))
    for pos in occ.keys():
        occ[pos]=collections.Counter()
    for line in listDic:
        occ[int(line['pos'])][line['type']]+=1
    mmOcc = pd.DataFrame(occ)
    return(mmOcc)


def make_occurrencies_tbl(seqDict):
    '''
    Make matrix of cooccurrencies of hexamer and primed sequence
    '''
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

# def makePWM(occTbl)

def makePWMtempl(counterTempl):
    '''
    From Counter of binding primers, makes PWM.
    '''
    count={0:collections.Counter(), 1:collections.Counter(), 2:collections.Counter(), 3:collections.Counter(), 4:collections.Counter(), 5:collections.Counter()}
    for k,v in counterTempl.items():
        for pos in count.keys():
            count[pos][k[pos]]+= v
    df = pd.DataFrame(count)
    df = df.drop([i for i in df.index if 'N' in i])
    return(df/df.sum())

def computeEntropy(pwm):
    '''
    Computes entropy of sequences binding with specific hexamer
    '''
    H=0
    for pos in pwm.iteritems():
        h=0
        for base in ["A","T","C", "G"]:
            if pos[1][base]!=0:
                h+= - pos[1][base]*math.log(pos[1][base], 4)
        H+=h
    return(H)


# filename='L1_R2_hex_entropy.txt'
# with open(filename, 'w') as f:
#     for row in tbl.iterrows():
#         print(row[0], computeEntropy(countBases(row[1])), sep='\t', file=f)
