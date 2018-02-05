import pysam as ps
import pandas as pd
import collections
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hexVSprimed import *
import pickle

bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/gk2a-2.sam'

with ps.AlignmentFile(bamfile,"rb") as bam:
    mmDic = mmPerQual(bam)
    # df=make_mm_table(mmDic)
    # df.to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mmPerPosition.csv')

with open('mmQual.pickle', 'wb') as handle:
    pickle.dump(mmDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
