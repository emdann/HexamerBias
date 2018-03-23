import pandas as pd
import pysam as ps
import numpy as np
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

def per_base_cov(params):
    '''
    Compute coverage for seq of interest based on template density
    '''
    seq, density = params
    posDic = {k:0 for k in range(len(seq))}
    for i in range(len(seq)-5):
        d = density[seq[i:i+6]]
        hexItems = [(k,v) for k,v in posDic.items() if k in range(i,i+6)]
        posDic.update(map(lambda kv: (kv[0], kv[1] + d), hexItems))
    return(pd.DataFrame(posDic, index=[0]).T)

def make_bed(bedEntry, fasta, threads=10):
    '''

    '''
    chr,start,end=bedEntry.split()
    seq = ps.FastaFile(fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
    workers = multiprocessing.Pool(threads)
    df = pd.DataFrame()
    for baseCov in workers.map(per_base_cov, [ (seq,density) for seq in [seq[i:i+100] for i in range(0,len(seq),100)]]):
            df = df.append(baseCov)
    df.index = range(int(start),int(end))
    covBed = []
    for row in df.iterrows():
        start,cov = row[0],row[1][0]
        covBed.append((chr, start, start+1, cov))
    return(pd.DataFrame(covBed, columns=['chr','start','end','coverage']))

# with ps.FastxFile(args.fasta) as chr:
#  	for entry in chr:
#  		seqs[entry.name] = entry.sequence.upper()

def running_mean(bed, win=100):
    '''
    Smoothen out base res artificial coverage
    '''
    bed.coverage = pd.rolling_mean(bed.coverage, window=win)
    bed = bed.dropna()
    return(bed)

def save_coverage_bed(bed, outfile):
    '''
    Save bed file of artificial coverage.
    '''
    with open(outfile, 'w') as out:
        print(bed.to_string(index=False, header=False), file=out)
    return('')


# abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz"
# covFile = "predictedCoverage_avgVAN1667.txt"
# abundance = pd.read_csv(abundanceFile, index_col=0, compression=findCompr(abundanceFile), header=None)
# coverage = pd.read_table(covFile, index_col=0, sep='\t',compression=findCompr(covFile))
#
# density = template_density(coverage.exp,abundance)
# bedEntry='chr6 52229197 52242854'
