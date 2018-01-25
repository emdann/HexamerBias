from hexVSprimed import *

mmfile = '/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_primed_seq.original.mismatch.txt'


newDic={}
with open(mmfile, 'rt') as f:
    f.readline()
    for line in f.readlines():
        newDic[line.split()[0]]=[line.split()[1], line.split()[2]]

kmerMM={}
for primedSeq,hex in newDic.values():
    if hex not in kmerMM.keys():
        kmerMM[hex]=[]
    kmerMM[hex].append(mismatchKmers((primedSeq,hex)))

outputfile='/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_avgMM.txt'
with open(outputfile,'w') as f:
    for k,i in kmerMM.items():
        print(k,np.mean(i), sep='\t', file=f)
