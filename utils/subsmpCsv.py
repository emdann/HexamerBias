import pandas
import random

def rand_subsample_csv(filename, sample_size = 20000, index_col=0):
    '''
    Take a random subsample of rows from big csv file'
    '''
    n = sum(1 for line in open(filename)) - 1 #number of records in file (excludes header)
    skip = sorted(random.sample(range(1,n+1),n-s)) #the 0-indexed header will not be included in the skip list
    df = pandas.read_csv(filename, skiprows=skip, index_col=index_col)
    return(df)
