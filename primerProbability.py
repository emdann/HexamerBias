import numpy as np
import pandas as pd
import itertools as it

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

def get_even_prob():
    seqs = [''.join(i) for i in list(it.product(list("ATCG"),repeat=6))]
    # Make structured array for concentrations
    x = np.array([[0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25],
                [0.25,0.25,0.25, 0.25, 0.25, 0.25]  ])
    # dtype=[('foo', 'i4'),('bar', 'f4'), ('baz', 'S10')])
    ppm = pd.DataFrame(x, index=["A", "T", "C", "G"])
    return(prob_from_ppm(ppm, seqs))

def get_proportional_coverage(dgTab, genomeAb):
    '''
    Compute predicted coverage with primer probabilities proportional to template abundance
    '''
    ppm = makePWMtempl(genomeAb)
    return(predictCoverage(dgTab, genomeAb, primer_ppm=ppm))
