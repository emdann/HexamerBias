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
            alPairs = r.get_aligned_pairs(with_seq=True)[:6]
            print(alPairs)
            for i in alPairs:
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

with open('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mmQualHex.pickle', 'wb') as handle:
    pickle.dump(mmDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('/Users/user/mnt/edann/hexamers/rnaseq/mmQualHex.pickle', 'rb') as handle:
    mmDic = pickle.load(handle)

counts = [len(i) for i in mmDic.values()]
widths = [c/float(sum(counts)) for c in counts]
plt.boxplot(mmDic.values(), labels=mmDic.keys())
plt.ylabel('Phread score')
plt.xticks(rotation=90)
plt.show()

plt.savefig('phred_score_hex.png')
