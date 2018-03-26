import pandas as pd
import pysam as ps
import numpy as np
import pyBigWig as pbw
import sys
from scipy import ndimage
from hexVSprimed import *
import multiprocessing


def template_density(covCol, abundance):
    '''
    Compute density (C/T) for each template sequence
    '''
    df  = pd.concat([abundance, covCol], axis=1, join_axes=[covCol.index])
    df.columns = ['abundance', 'coverage']
    density = df.coverage/df.abundance
    return(density)

def splitDataFrameIntoSmaller(df, chunkSize = 10000):
    listOfDf = list()
    numberChunks = len(df) // chunkSize + 1
    for i in range(numberChunks):
        listOfDf.append(df[i*chunkSize:(i+1)*chunkSize])
    return listOfDf

def per_base_cov(params):
    '''
    Compute coverage for seq of interest based on template density
    '''
    seq, density, start = params
    posDic = {k:np.float64() for k in range(start, start+len(seq))} # Specify to have same type in the end
    for i in range(len(seq)-5):
        if 'N' not in seq[i:i+6]: ## No density on sites with N but the position is there
            d = density[seq[i:i+6]]
            hexItems = [(k,v) for k,v in posDic.items() if k in range(start+i,start+i+6)]
            posDic.update(map(lambda kv: (kv[0], kv[1] + d), hexItems))
    return(posDic)

def running_mean_dic(baseCov, win=100):
    '''
    Smoothen out base resolution artificial coverage with running average
    '''
    bed = pd.DataFrame(baseCov, index=[0]).T
    bed.columns=['coverage']
    bed.coverage = round(bed.coverage.rolling(window=win).mean(),4) # reducing significant digits to be able to collapse
    bed = bed.dropna()
    return(bed.to_dict()['coverage'])

def kernel_smoothing(baseCov, sigma=20):
    '''
    gaussian kernel smoothing
    '''
    bed = pd.DataFrame(baseCov, index=[0]).T
    bed.columns=['coverage']
    bed.coverage = ndimage.gaussian_filter1d(bed.coverage, sigma=sigma)
    return(bed.to_dict()['coverage'])


### Functions to make BigWig ###

def make_BigWig_header(refgen):
    '''
    Make header for bigWig file with chromosome length from fasta of refgen
    '''
    header=[]
    with ps.FastxFile(refgen) as f:
        for entry in f:
            header.append((entry.name, len(entry.sequence)))
    def get_chr_num(chr):
        return(chr.split('chr')[1])
    return(sorted(header, key=lambda x: get_chr_num(x[0])))

def add_seq_to_bigWig(seq,chrom,start,density, bw_with_header, filename, smoothFunction = 'kernel', threads=10):
    '''
    '''
    workers = multiprocessing.Pool(threads)
    for posDic in workers.map(per_base_cov, [ (seq,density,start+s) for seq,s in [(seq[i:i+1099],i) for i in range(0,len(seq),1000)]]):
            if smoothFunction=='running avg':
                posDic = running_mean_dic(posDic, win=100)
            if smoothFunction=='kernel':
                posDic = kernel_smoothing(posDic, sigma=20)
            else:
                print('Unrecognized smoothing function! Please choose "running avg" or "kernel"')
                return ''
            bw_with_header.addEntries(chrom, [k for k in posDic.keys()], values=[v for v in posDic.values()], span=1)
    return(bw_with_header)

def save_bigWig(beds,refgen_fasta,density, outfile, threads=10):
    '''
    Input: list of bed entries (str of one line), fasta of reference genme, genome file with chrom sizes
    Output:
    '''
    bw = pbw.open(outfile, 'w')
    bw.addHeader(make_BigWig_header(refgen_fasta))
    for entry in beds:
        chr,start,end = entry.split()
        print('Processing entry ', entry)
        seq = ps.FastaFile(refgen_fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
        add_seq_to_bigWig(seq, chr,int(start), density, bw, outfile, threads=threads)
    bw.close()
    return(bw)

### Functions to make BedGraph ###

def make_bed(bedEntry, fasta, density, threads=10):
    '''

    '''
    chr,start,end=bedEntry.split()
    seq = ps.FastaFile(fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
    workers = multiprocessing.Pool(threads)
    df = pd.DataFrame()
    for baseCov in workers.map(per_base_cov, [ (seq,density) for seq in [seq[i:i+100] for i in range(0,len(seq),100)]]):
            df = df.append(pd.DataFrame(baseCov, index=[0]).T)
    df.index = range(int(start),int(end))
    covBed = []
    for row in df.iterrows():
        start,cov = row[0],row[1][0]
        covBed.append((chr, start, start+1, cov))
    return(pd.DataFrame(covBed, columns=['chr','start','end','coverage']))

def collapse_coverage_bed(bed):
    '''
    Collapse base resolution bed file for coverage if the coverage field
    of consecutive entries is the same
    '''
    c,s,e,d = bed.iloc[0]
    collapsedBed = pd.DataFrame()
    for entry in bed.iterrows():
        chr,start,end,coverage = entry[1]
        if coverage == d:
            c,s,e,d = c,s,end,d
        else:
            df = pd.DataFrame({'chr':c,'start':s,'end':e,'coverage':d}, index=[0])
            collapsedBed = collapsedBed.append(df)
            c,s,e,d = entry[1]
    df = pd.DataFrame({'chr':c,'start':s,'end':e,'coverage':d}, index=[0])
    collapsedBed = collapsedBed.append(df)
    collapsedBed = collapsedBed[['chr','start','end','coverage']]
    return(collapsedBed)

def running_mean_bed(bed, win=100, threads=10):
    '''
    Smoothen out base res artificial coverage
    '''
    bed.coverage = round(bed.coverage.rolling(window=win).mean(),4) # reducing significant digits to be able to collapse
    bed = bed.dropna()
    listBed = splitDataFrameIntoSmaller(bed, chunkSize=100)
    workers = multiprocessing.Pool(threads)
    collapsedBed = pd.DataFrame()
    for mergedBed in workers.map(collapse_coverage_bed, [ bed for bed in listBed]):
            collapsedBed = collapsedBed.append(mergedBed)
    return(mergedBed)

def artificial_cov(beds,fasta,density, threads=10, win=100):
    '''
    Computes artificial coverage with running mean on bed entries
    '''
    out_bed=pd.DataFrame()
    for entry in beds:
        print("Processind entry " + entry)
        covBed = make_bed(entry,fasta,density, threads = threads)
        if win!=0:
            print("Running average on entry " + entry)
            covBed = running_mean(covBed, win=win, threads=threads)
        out_bed = out_bed.append(covBed)
    return(out_bed)

def save_coverage_bed(bed, outfile):
    '''
    Save bed file of artificial coverage.
    '''
    with open(outfile, 'w') as out:
        print(bed.to_string(index=False, header=False), file=out)
    return('')

# abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz"
# covFile = "predictedCoverage_avgVAN1667.txt"
