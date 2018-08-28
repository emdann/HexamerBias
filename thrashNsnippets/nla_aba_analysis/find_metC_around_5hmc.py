import pysam as ps
import fnmatch
import collections
import pandas as pd

cOHfile='CG-abaBS-1000ES-3n1_lmerged_R1.CTOB.uniq5hmc.bed'
hydroxySites=[]
with open(cOHfile, 'r') as f:
    for line in f.readlines():
        hydroxySites.append(line)

bed.flank(g='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.genome', b=1000)
