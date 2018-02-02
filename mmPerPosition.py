import pysam as ps
import pandas as pd
import collections
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hexVSprimed import *

bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/gk2a-2.sam'

with ps.AlignmentFile(bamfile,"rb") as bam:
    hexmm,countedHex = mmPerHex(bam),countHex(bam)
    df = pd.DataFrame([hexmm, countedHex]).T
    df.columns=['avgMM', 'usage']
    df.to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/hexMMusage.csv')
