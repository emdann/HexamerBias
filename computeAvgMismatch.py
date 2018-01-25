from hexVSprimed import *
import collections


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
    kmerMM[hex].append(mismatchKmers((primedSeq,hex), tolerateBSmm=True))

outputfile='/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_avgMM_tolerateBSmm.txt'
with open(outputfile,'w') as f:
    for k,i in kmerMM.items():
        print(k,np.mean(i), sep='\t', file=f)

# hexs=[]
# primedregs=[]
# with open(mmfile, 'rt') as f:
#     f.readline()
#     for line in f.readlines():
#         hexs.append(line.split()[2])
#         primedregs.append(line.split()[1])
#
# hexCount = collections.Counter(hexs)
# for k in list(hexCount.keys()):
#   if 'N' in k:
#     hexCount.pop(k)
#
# primedCount = collections.Counter(primedregs)
# for k in list(primedCount.keys()):
#   if 'N' in k:
#     primedCount.pop(k)
#
#
# file='/hpc/hub_oudenaarden/edann/hexamers/mismatch/L1_R1_hexVSprimed.count.txt'
# with open(file, 'w') as f:
#     print('seq', 'primed_region_count', 'usage_count', sep='\t',file=f)
#     for el in primedCount.keys():
#         print(el, primedCount[el], hexCount[el], sep='\t', file=f)
