import pysam as ps
import fnmatch
import collections
import pandas as pd

# fasta='/hpc/hub_oudenaarden/cgeisenberger/data_analysis/reBS/CG-abaBS-1000xES-3n1/CG-abaBS-1000ES-3n1_AHFKWYBGX5_S13_L001_R1_001.fastq.gz'

def strip_barcode(readSeq, lastBarcode):
    '''
    Removes the adapter barcode part of the fastq reads
    '''
    # print(readSeq.find(lastBarcode))
    # strippedRead = readSeq[readSeq.find(lastBarcode) + len(lastBarcode):]
    bcPos = readSeq.find(lastBarcode)
    if bcPos != -1:
        strippedRead = readSeq[bcPos+len(lastBarcode):]
        return(strippedRead)
    else:
        return

def base_comp_after_cutsite(seqs):
    '''
    Find distance between C-OH and cutsite in abaSeq
    '''
    posBaseComp = {}
    length = max([len(s) for s in seqs])
    for i in range(length):
        p=[s[i:i+1] for s in seqs if i in range(len(s))]
        posBaseComp[i] = collections.Counter(p)
    baseTab = pd.DataFrame(posBaseComp)
    return(baseTab)

def make_basecomp_table_aba(fastqfile):
    seqs=[]
    lastBarcode = "AACCAGAACAGTGGTATCAACGCAGAGTAACCAACG"
    with ps.FastxFile(fastqfile) as f:
        for entry in f:
            seqs.append(strip_barcode(entry.sequence, lastBarcode))
    seqs = [i for i in seqs if i is not None]
    outfile = fastqfile.split('/')[-1].split('.fastq')[0] + '.basematrix.csv'
    df = base_comp_after_cutsite(seqs)
    df.to_csv(outfile)
    return(df)

def make_basecomp_table_aba_fromfile(file):
    seqs=[]
    with open(file, 'r') as f:
        for entry in f.readlines():
            seqs.append(entry.rstrip())
    seqs = [i for i in seqs if i is not None]
    outfile = file.split('/')[-1].split('.txt')[0] + '.basematrix.csv'
    df = base_comp_after_cutsite(seqs)
    df.to_csv(outfile)
    return(df)

def metcall_table_aba(file):
    seqs163=[]
    seqs147=[]
    with open(file, 'r') as f:
        for entry in f.readlines():
            entry=entry.strip()
            flag,seq,metCall = entry.split(' \t ')
            if flag=='163':
                seqs163.append(metCall)
            if flag=='147':
                seqs147.append(metCall)
    outfile147 = file.split('/')[-1].split('.txt')[0] + '.CTOT.basematrix.csv'
    outfile163 = file.split('/')[-1].split('.txt')[0] + '.CTOB.basematrix.csv'
    df147 = base_comp_after_cutsite(seqs147)
    df163 = base_comp_after_cutsite(seqs163)
    df147.to_csv(outfile147)
    df163.to_csv(outfile163)
    return(df163,df147)
