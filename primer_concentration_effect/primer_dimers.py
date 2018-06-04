import pysam as ps
import itertools as it
import collections

reads='CG-pbat-gDNA-zfishTotal-noBS-1preAmp-handMix_lmerged_R1.fastq.gz'

lenCount = collections.Counter()
with ps.FastxFile(reads) as fastq:
    for entry in fastq:
        lenCount[len(entry.sequence.upper())]+=1
