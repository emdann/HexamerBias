import pysam as ps
import pandas as pd
import collections
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hexVSprimed import *
import pickle

bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/zebrafish/gk2a-2.sam'

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

with ps.AlignmentFile(bamfile,"rb") as bam:
    mmDic = mmPerQual(bam)
    # df=make_mm_table(mmDic)
    # df.to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mmPerPosition.csv')

with open('mmQual.pickle', 'wb') as handle:
    pickle.dump(mmDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
