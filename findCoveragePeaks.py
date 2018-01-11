import pandas as pd
import peakutils
import numpy as np
from matplotlib import pyplot
import multiprocessing


file='/Users/user/mnt/edann/VAN1667_depth.bed'


def findPeaks(cov):
    """
    detects peaks in depth.bed file (use small chunks)
    returns df in bed format of 1 kb region around detected peaks
    """
    chr=np.array([])
    x=np.array([])
    y=np.array([])
    chrom=cov.chrom[0]
    for i in range(min(cov.start), max(cov.end)):
        df=cov[(cov.start <= i) & (cov.end > i)]
        if df.empty:
            x=np.append(x, i)
            y=np.append(y, 0)
        else:
            x=np.append(x,i)
            y=np.append(y,df.depth)
            chrom=df.chrom
        chr=np.append(chr, chrom)
    complCov=pd.DataFrame({"pos":x,"depth":y, "chrom":chr})
    complCov['runmeanDepth']=pd.rolling_mean(complCov.depth, 500)
    pyplot.plot(complCov.pos, complCov.runmeanDepth)
    pyplot.plot(complCov.pos, complCov.depth)
    pyplot.show()
    indexes = peakutils.indexes(complCov.runmeanDepth, thres=0.5, min_dist=1000)
#    pplot(complCov.pos, complCov.runmeanDepth, indexes)
#    pyplot.show()
    peaks = complCov.iloc()[indexes]
    peaksBed = pd.DataFrame({'chrom':peaks.chrom, 'start':peaks.pos-500, 'end':peaks.pos+500, 'cov':peaks.depth})
    return(peaksBed)

def read_in_chunks(file, chunk_size=1000):
    """Lazy function (generator) to read a file piece by piece.
        Default chunk size: 1k."""
    firstrow = 0
    while True:
        data = pd.read_csv(file, header=None, skiprows=firstrow, nrows=chunk_size, names=['chrom', 'start', 'end', 'depth' ], sep='\t')
        firstrow = firstrow+chunk_size
        if data.empty:
            break
        yield data


workers = multiprocessing.Pool(10)
p=pd.DataFrame([])
for peaks in workers.map(findPeaks, [ cov for cov in read_in_chunks(file)]):
    p=p.append(peaks)

p.to_csv('VAN1667_depth_peaks.csv')
