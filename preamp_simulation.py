import pandas as pd
import pysam as ps
import numpy as np
import random
import sys
import multiprocessing
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

def seq_from_bed_entry(params):
    bedEntry,refgen_fasta=params
    chr,start,end = bedEntry.split()
    seq = ps.FastaFile(refgen_fasta).fetch(reference=chr, start=int(start), end=int(end)).upper()
    return(seq)

def fragment_entry(params):
    '''
    Takes bed entry and fragments the corresponding regions (as BS conversion)
    Output: dict of {coords:seq}
    '''
    bedEntry,refgen_fasta=params
    seq = seq_from_bed_entry((bedEntry, refgen_fasta))
    chr,start,end = bedEntry.split()
    p = 0
    fragmentPool = {}
    while p < len(seq):
        pEnd=p+int(np.random.normal(loc=500, scale=200))
        fragment = seq[p:pEnd]
        frag_start = int(start) + p
        frag_end = int(end) + pEnd
        coords = ":".join([chr,str(frag_start),str(frag_end)])
        fragmentPool[coords+":+"] = fragment
        fragmentPool[coords+":-"] = str(Seq(fragment, generic_dna).reverse_complement())
        p=pEnd
    fragmentPoolSizeFilt = {coord:f for coord,f in fragmentPool.items() if len(f) > 200}    # Exclude fragments under set length
    return(fragmentPoolSizeFilt)

def find_kmers(params):
    string, k, addTo = params
    string = string.upper()
    kmerCounts = collections.Counter() if addTo is None else addTo
    n = len(string)
    for i in range(0, n-k+1):
        kmerCounts[ string[i:i+k] ] += 1
    return kmerCounts

def abundance_from_fragment_pool(kmerFragment):
    '''
    Input: dict of collection.Counters for kmer abundance in each fragment of the genome
    Output: dataframe of total template abundance
    '''
    abund = pd.DataFrame(kmerFragment).T.sum()
    abundDf = pd.DataFrame(abund)
    abundDf.index.name='template'
    abundDf.reset_index(inplace=True)
    abundDf = abundDf.rename(index=str, columns={abundDf.columns[1]:'abundance'})
    return(abundDf)

def predict_coverage(T, primerProb, keqs, eps, nreads=1000000):
    keqsDf = keqs.merge(primerProb, on='primer').merge(T, on='template')
    keqsDfPhi = keqsDf.assign(phi=keqsDf.primerProb*keqsDf.keq).groupby('template').agg({'phi':'sum', 'abundance':'max'})
    predCovDf = keqsDfPhi.assign(predCov=keqsDfPhi.abundance*eps*(keqsDfPhi.phi/(1+keqsDfPhi.phi)))
    normPredCov =  predCovDf.assign(normPredCov=predCovDf.predCov*nreads/sum(predCovDf.predCov))
    dens=predCovDf.assign(dens=normPredCov.normPredCov/normPredCov.abundance)
    dens.reset_index(inplace=True)
    return(dens)

def fragment_amplification(params):
    fragment, dens,coord = params
    fragmentDf = pd.DataFrame(fragment, index=['kmerAb']).T
    fragmentDf.index.name='template'
    fragmentDf.reset_index(inplace=True)
    fragDens = fragmentDf.dropna().merge(dens, how='left', on='template')
    nAmp=0
    for temp in fragDens.iterrows():
        template=temp[1]
        n=0
        while n<template.kmerAb:
            rand = random.uniform(0,1)
            if rand<template.dens:
                nAmp +=1
            n += 1
    return(coord,nAmp)

def multiply_kmer_fragments(fragCounter, kmerFragment):
    kmerFragmentsNew = {}
    for fr,count in fragCounter.items():
        newKmer = {k:(v*count) for k,v in kmerFragment[fr].items()}
        kmerFragmentsNew[fr]=newKmer
    return(kmerFragmentsNew)

def preamp_round(kmerFragment, fragCounter, primerProb, keqs, eps, nreads=1000000, threads=5):
    ## Predict coverage
    dens = predict_coverage(abundance_from_fragment_pool(kmerFragment), primerProb, keqs, eps, nreads=nreads)
    ## Amplify
    workers = multiprocessing.Pool(threads)
    amplifiedCounter = fragCounter.copy()
    for (coord,nAmp) in workers.imap_unordered(fragment_amplification, [[fr,dens,coord] for coord,fr in kmerFragment.items()]):
        amplifiedCounter[coord]+=nAmp
    # for fr in fragCounter.keys(): # <---- PARALLELIZABLE
    #     nAmp=fragment_amplification(kmerFragment[fr], dens)
    #     amplifiedCounter[fr]+=nAmp
    return(amplifiedCounter)

def save_preamp_result(amp, outfile):
    #     header = make_BigWig_header(refgen_fasta)
    # bw = pbw.open(outfile, 'w')
    # header=make_BigWig_header(refgen_fasta)
    # chroms = [e[0] for e in header]
    # bw.addHeader(header)
    # print('--- Saving BigWig ---', flush=True)
    # for chrom in chroms:
    #     print('Writing entries from '+ chrom, flush=True)
    #     exChr = sorted([i for i in amp.keys() if i[0]==chrom])
    #     if len(exChr)>0:
    #         bw.addEntries(chrom, [k[1] for k in exChr], values=[k[2] for k in exChr], span=1)
    # bw.close()
    # outputDf = pd.DataFrame(amp1, index=['n']).T.reset_index().rename(index=str, columns={'index':'reg'})
    # outputDf.assign(chr=[reg.split(':')[0] for reg in outputDf.reg], start=[reg.split(':')[1] for reg in outputDf.reg], end=[reg.split(':')[2] for reg in outputDf.reg])
    with open(outfile, 'w') as out:
        for coord,n in amp1.items():
            chr,start,end,strand = coord.split(':')
            print(chr,start,end,".", n, strand, sep='\t', file=out)
    return(None)

# ---------------------------------------------------------- #

refgen_fasta='/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'
bedfile='/hpc/hub_oudenaarden/edann/preAmp_rounds_VAN2590/random_10kb.42.mm10.bed'
keqsFile='/hpc/hub_oudenaarden/edann/WGS_keqs.csv'
threads=5

bedList=[]
with open(bedfile, 'r') as f:
    for entry in f.readlines():
        bedList.append(entry)


# genome={}
# with ps.FastxFile(genomefile) as f:
#     for entry in f:
#         genome[entry.name]=entry.sequence

# randGenome = entry.sequence[10000000:10500000].upper()

## Model parameters
primerProb = get_even_prob()
primerProb.index.name='primer'
primerProb.reset_index(inplace=True)
# T1 = find_kmers((randGenome, 6))    # Count initial total kmer abundance
# # ADD STRAND SPECIFICITY!
# T1df = pd.DataFrame(T1, index=['abundance']).T
# T1df.index.name='template'
# T1df.reset_index(inplace=True)
keqs = pd.read_csv(keqsFile)
eps = 500

## Make random breaks
workers = multiprocessing.Pool(threads)
fragmentPoolFull={}
for entryPool in workers.imap_unordered(fragment_entry, [(bed, refgen_fasta) for bed in bedList]):
    fragmentPoolFull.update(entryPool)



## Count kmers in every fragment
n = 0
kmerFragment = {}
for coord,frag in fragmentPoolFull.items():
    kmers = find_kmers((frag, 6, None))
    kmerFragment[coord] = kmers

fragCounter = collections.Counter(fragmentPoolFull.keys())

# smpFrags=random.sample(list(kmerFragment),10)
# smpFragsKmers = {x:kmerFragment[x] for x in smpFrags}
# smpFragsCounter = {x:fragCounter[x] for x in smpFrags}

## Amplification rounds
amp1 = preamp_round(kmerFragment, fragCounter, primerProb, keqs, eps, nreads=1000000)
amp2 = preamp_round(multiply_kmer_fragments(amp1, kmerFragment), amp1, primerProb, keqs, eps, nreads=1000000)
amp3 = preamp_round(multiply_kmer_fragments(amp2, kmerFragment), amp2, primerProb, keqs, eps, nreads=1000000)
amp4 = preamp_round(multiply_kmer_fragments(amp3, kmerFragment), amp3, primerProb, keqs, eps, nreads=1000000)
amp5 = preamp_round(multiply_kmer_fragments(amp4, kmerFragment), amp4, primerProb, keqs, eps, nreads=1000000)

save_preamp_result(amp1, 'simulation_preAmp1.bed')
save_preamp_result(amp2, 'simulation_preAmp2.bed')
save_preamp_result(amp3, 'simulation_preAmp3.bed')
save_preamp_result(amp4, 'simulation_preAmp4.bed')
save_preamp_result(amp5, 'simulation_preAmp5.bed')
