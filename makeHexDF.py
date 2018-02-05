import pysam as ps
import pandas as pd
import collections
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hexVSprimed import *

bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/zebrafish/gk2a-2.sam'
fasta='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2_primed_seq.fa'

templDic={}
with ps.FastxFile(fasta) as templ:
	for entry in templ:
		seq, name = entry.sequence.upper(), entry.name
		templDic[name]=seq

with ps.AlignmentFile(bamfile,"rb") as bam:
    UsageMmDicPrimer = usageMmPerHex(bam, ref='primer')
    UsageMmDicTempl = usageMmPerHex(bam, ref='template', templDic=templDic)

dfPrimer = pd.DataFrame(UsageMmDicPrimer).T
dfPrimer.columns = ['bindingEventsPrimer','mmEventsPrimer', 'avgMmPrimer']

dfTempl = pd.DataFrame(UsageMmDicTempl).T
dfTempl.columns = ['bindingEventsTempl','mmEventsTempl', 'avgMmTempl']

df = pd.concat([dfTempl, dfPrimer], axis=1)
df.to_csv('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2.hexTable.csv')
