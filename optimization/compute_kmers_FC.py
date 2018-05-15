#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Compute FC between kmers abundance of region of interest and random regions in genome. By Emma Dann")
argparser.add_argument('kmersROI', type=str, help='path to file of kmer abundance in regions of interest')
argparser.add_argument('kmersRandom', type=str, help='path to file of kmer abundance in random regions')
args = argparser.parse_args()

def compute_kmers_FC(kmerFile, randomKmerFile, save_csv=True):
    '''
    Compute log2(FC) between kmer count in regions of interest and random regions
    to be used for DE optimization
    Input:
        kmerFile: path to file of kmer abundance in regions of interest
        randomKmerFile: path to file of kmer abundance in random regions
        save_csv: logical for wheather the csv of kmers should be saved in same folder as kmerFile
    '''
    targetKmers = pd.read_csv(kmerFile, index_col=0, header=None)
    randomKmers = pd.read_csv(randomKmerFile, index_col=0, header=None)
    enrichKmers = pd.concat([randomKmers, targetKmers], axis=1)
    enrichKmers.columns = ['random','target']
    foldChange = np.log2(enrichKmers.target/enrichKmers.random).dropna()
    if save_csv:
        foldChange.to_csv(kmerFile.rstrip('.kmersTot.csv')+'.kmersFC.csv')
    return(foldChange)

compute_kmers_FC(args.kmersROI, args.kmersRandom, save_csv=True)
