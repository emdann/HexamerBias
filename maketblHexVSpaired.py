from hexVSprimed import *
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Get hexamers used in fasta file.\n By Emma Dann")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('primedreg', type=str, help='Fasta input')
args = argparser.parse_args()

hexFile=args.fasta
primedregFile=args.primedreg

dic={}
with ps.FastxFile(primedregFile) as fastq:
    for entry in fastq:
        dic[entry.name.split('_')[0]]=[entry.sequence.upper()]

if hexFile.endswith(('.fasta','.fastq','.fa', '.fa.gz', '.fasta.gz','.fastq.gz')):
    with ps.FastxFile(hexFile) as fastq:
        for entry in fastq:
    #        print(entry.name)
            if entry.name in dic.keys():
                dic[entry.name].append(entry.sequence[:6])
else:
    with open(hexFile, 'r') as f:
        for line in f.readlines():
            l=line.strip().split("\t")
            dic[l[0]].append(l[1])

with open(primedregFile.split("/")[-1].split('.fa')[0] +'.mismatch.txt', 'w') as output:
   print('read', 'primedSeq', 'hex', sep='\t', file=output)
   for key,val in dic.items():
       print(key, val[0], val[1], sep='\t', file=output)

# newDic={}
# with open('/hpc/hub_oudenaarden/edann/hexamers/rnaseq/'+primedregFile.split("/")[-1].split('.fa')[0] +'.mismatch.txt', 'rt') as f:
#     f.readline()
#     for line in f.readlines():
#         newDic[line.split()[0]]=[line.split()[1], line.split()[2]]

tbl = make_occurrencies_tbl(dic)
tbl.to_csv(primedregFile.split("/")[-1].split('.fa')[0] +'.csv')
