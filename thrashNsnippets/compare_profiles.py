import argparse
import sys
sys.path.insert(0,'/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage')
from cov_from_density import *



## Comparing two profiles from bigWig files
a = pbw.open('mm10.random.40.artCov.realProportional.bw')
b = pbw.open('mm10.random.40.artCov.pooled.kernel.bw')

dic1 = dict(((el[0],el[2]) for el in b.intervals('chr1')))
dic2 = dict(((el[0],el[2]) for el in a.intervals('chr1')))

sharedKeys = dic1.keys() & dic2.keys()

## Better with arrays
intervals = []
for chrom in a.chroms().keys():
    ints = a.intervals(chrom)
    if ints is not None:
        intervals.extend(ints)
arr1 = np.array(intervals)
arr2 = np.array(b.intervals('chr2'))
rho,pval = scipy.stats.spearmanr(arr1.T[2], arr2.T[2])
