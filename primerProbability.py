import numpy as np
import pandas as pd
'''
Compute the probability for all possible hexamers base on manufacturers concentrations.
Input: matrix of base concentration per position
Output: probability for all hexamers
'''

sequences = list(pd.read_csv(abfile, index_col=0).index)


# Make structured array for concentrations
x = np.array([[0.5,0.25,0.25, 0.25, 0.25, 0.25],
            [0.0,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25],
            [0.25,0.25,0.25, 0.25, 0.25, 0.25]  ])
# dtype=[('foo', 'i4'),('bar', 'f4'), ('baz', 'S10')])
ppm = pd.DataFrame(x, index=["A", "T", "C", "G"])

def prob_from_ppm(ppm, sequences):
    '''
    Compute probability of each sequence given a positional probability matrix
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
