import numpy as np
import pandas as pd
import itertools as it
from hexVSprimed import *
import random

def prob_from_ppm(ppm, sequences):
    '''
    Compute the probability for all possible hexamers based on manufacturers concentrations.
    Input: matrix of base concentration per position
    Output: probability for all hexamers
    '''
    seqProb = {}
    for seq in sequences:
        prob = 1
        for i in range(len(seq)):
            prob *= ppm.loc[seq[i],i]
        seqProb[seq] = prob
    seqProbDf = pd.DataFrame(seqProb, index=[0]).T
    seqProbDf.columns=['primerProb']
    return(seqProbDf.sort_values(by=["primerProb"]))

def from_vec_to_ppm(vector):
    '''
    Transform np.array line of probabilities (every 4 we have one position) in ppm
    '''
    arr = vector.reshape(4,int(vector.size/4))
    ppm = pd.DataFrame(arr, index=["A", "T", "C", "G"])
    return(ppm)

def all_hexamers():
    seqs = [''.join(i) for i in list(it.product(list("ATCG"),repeat=6))]
    return(seqs)

def get_even_prob():
    # Make structured array for concentrations
    x = np.array([[0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25]  ])
    # dtype=[('foo', 'i4'),('bar', 'f4'), ('baz', 'S10')])
    ppm = pd.DataFrame(x, index=["A", "T", "C", "G"])
    seqs = all_hexamers()
    return(prob_from_ppm(ppm, seqs))

def get_proportional_coverage_ppm(dgTab, genomeAb):
    '''
    Compute predicted coverage with primer probabilities proportional to template abundance
    '''
    ppm = makePWMtempl(genomeAb)
    return(predictCoverage(dgTab, genomeAb, primer_ppm=ppm))

def get_proportional_probability(genomeAb):
    '''
    Compute predicted coverage with primer probabilities proportional to template abundance
    '''
    proportionalProb = genomeAb/genomeAb.sum()
    probs = pd.DataFrame(proportionalProb)
    probs.columns=['primerProb']
    return(probs)

def heaviside_step_function(r,r0=0.50, m=1, nn=5):
    if m>0 and m<nn:
        r=r
    if m>nn and m<3*nn:
        r=(r/2)*(1/r)
    if m>3*nn:
        r = ((r/2)*(1/r))* np.sin(2*r0*3*nn)
    return(r)

def random_nucleotide_probability(step):
    '''
    '''
    probs = []
    totalProb = 0
    # Choose a random probability p1 from a uniform distribution in
    # the range (0, 1), then choose p2 in the range (0, 1 - p1), etc.
    for i in range(3):
        up_to = 1 - totalProb
        p = random.choice(np.arange(0.0 , up_to+(step*0.01) ,step))  # putting +0.0001 to include the 1
        probs.append(p)
        totalProb += p
    probs.append(1 - totalProb)
    bases = ['A', 'C', 'G', 'T']
    random.shuffle(bases)
    new_prob_dict = {}
    for base, prob in zip(bases, probs):
        new_prob_dict[base] = prob
    return(pd.DataFrame(new_prob_dict, index=[0]).T)

def random_ppm(step, n=6):
    '''
    '''
    ppm=pd.DataFrame([], index=["A", "T", "C", "G"])
    for pos in range(n):
        ppm = pd.concat([ppm, random_nucleotide_probability(step=step)], axis=1)
    ppm.columns=range(n)
    return(ppm)

def change_nucleotide_probability(prob_series, base_to_optimize, step=0.25):
    '''
    '''
    new_prob = pd.Series(index=["A", "T", "C", "G"]) # avoiding renaming and not reassigning
    new_prob[base_to_optimize] = prob_series[base_to_optimize] + step
    mask = new_prob.index.isin([base_to_optimize])
    totalProb = new_prob[base_to_optimize]
    probs = []
    for i in range(2):
        up_to = 1 - totalProb
        p = random.choice(np.arange(0.0 , up_to+(abs(step)*0.01) ,abs(step)))  # putting +0.0001 to include the 1
        probs.append(p)
        totalProb += p
    probs.append(1 - totalProb)
    bases = [base for base in new_prob.index if base is not base_to_optimize]
    random.shuffle(bases)
    for base, prob in zip(bases, probs):
        new_prob[base] = prob
    return(new_prob)
