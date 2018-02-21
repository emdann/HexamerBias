import pandas as pd
import pysam as ps

def num_reads_per_cell(bamfile, to_file = False):
    '''
    Makes dictionary of number of reads for each cell
    '''
    
