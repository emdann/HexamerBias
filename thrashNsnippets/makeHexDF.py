import pysam as ps
import pandas as pd
import collections
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hexVSprimed import *
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Extract position of primer placement from trimmed section of aligned reads. By Emma Dann")
argparser.add_argument('bam', type=str, help='Input bam file')
argparser.add_argument('primedseq_fasta', type=str, help='Fasta file of primed regions')
argparser.add_argument('abundance', type=str, help='Hex abundance in refgenome')
args = argparser.parse_args()

bamfile = args.bam
fasta = args.primedseq_fasta
abFile =  args.abundance

# bamfile='/hpc/hub_oudenaarden/aalemany/emma-adi/mouse/SvdB11d4-MitoTrackerThird-Satellites-Adult.sam'
# fasta='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/SvdB11d4-MitoTrackerThird-Satellites-Adult.sam_primed_reg.qname.bed.fa'

templDic={}
with ps.FastxFile(fasta) as templ:
	for entry in templ:
		seq, name = entry.sequence.upper(), entry.name
		templDic[name]=seq

with ps.AlignmentFile(bamfile,"rb") as bam:
    UsageMmDicPrimer = usageMmPerHex(bam, ref='primer')

with ps.AlignmentFile(bamfile,"rb") as bam:
    UsageMmDicTempl = usageMmPerHex(bam, ref='template', templDic=templDic)

dfPrimer = pd.DataFrame(UsageMmDicPrimer).T
dfPrimer.columns = ['bindingEventsPrimer','mmEventsPrimer', 'avgMmPrimer']

dfTempl = pd.DataFrame(UsageMmDicTempl).T
dfTempl.columns = ['bindingEventsTempl','mmEventsTempl', 'avgMmTempl']

df = pd.concat([dfTempl, dfPrimer], axis=1)

# Add abundance
abDf = pd.read_csv(abFile, sep='\t')
abDf.columns = ['abundance', 'kmer']
abDf.index = abDf.kmer
df2 = pd.concat([df,abDf.abundance], axis=1, join='inner')
df3 = df2.drop([i for i in df2.index if 'N' in i])

# Save
output = '/hpc/hub_oudenaarden/edann/hexamers/rnaseq/'+bamfile.split('/')[-1][:-4]+'.phredFilter.hexTab.wAb.csv'
df3.to_csv(output)
