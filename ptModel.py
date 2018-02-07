from hexVSprimed import *
import pysam as ps
import collections
import argparse
import numpy as np
import pandas as pd
import multiprocessing

def extractDeltaG(complex, primer, template):
    '''
    Extract predicted DeltaG for specific primer-template interaction
    based on concentration of primer-template complex, of primer in the experiment,
    of number of template sequences in reference genome.
    '''
    kb = ...
    T = ...
    dG = np.log(complex/(primer*template)) * (kb*T)
    return(dG)
