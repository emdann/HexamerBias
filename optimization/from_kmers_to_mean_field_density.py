import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization')
from primerProbability import *
from predictCovBs import *
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *

### WHAT AM I TRYING TO DO?
# Given a set of regions for which we want to maximize coverage, we can compute kmer
# abundance and multiply the abundances by the density. This is what we want to maximize

def mean_field_density(kmerCounts, density):
    '''
    From kmer abundance of regions we want to maximize, compute total coverage density
    multiplying kmer abundance by density.
    '''
    df = pd.concat([kmerCounts, density], axis=1)
    mul=df[1]*df[0]
    return(mul.sum())

kmerCountsFiles = "CTCF.flank60.kmers.csv"
kmerCounts = pd.read_csv(kmerCountsFiles, index_col=0, header=None)
