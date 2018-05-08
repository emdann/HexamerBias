#!/usr/bin/env python3
import pysam
import argparse
import collections
import multiprocessing
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Count kmers in fasta file. By Buys de Barbanson")
argparser.add_argument('fasta', type=str, help='Fasta input')
argparser.add_argument('-s', type=str, default='+', required=False, help='strand')
argparser.add_argument('-k', type=int, default=6, required=False, help='Kmer size')
argparser.add_argument('-t', type=int, default=8, required=False, help='Amount of threads to use.')
args = argparser.parse_args()

def find_kmers(params):
    string, k, addTo = params
    string = string.upper()
    kmerCounts = collections.Counter() if addTo is None else addTo
    n = len(string)
    for i in range(0, n-k+1):
        kmerCounts[ string[i:i+k] ] += 1
    return kmerCounts

strand=args.s
threads=args.t
k=args.k
fasta=args.fasta

workers = multiprocessing.Pool(threads)
finalKmerCounts = collections.Counter()
if strand=='+':
    with pysam.FastxFile(fasta) as f:
        for kmerCounts in workers.imap_unordered(find_kmers, [ (str(entry.sequence), k, None) for entry in f]):
            finalKmerCounts+=kmerCounts
elif strand=='-':
    with pysam.FastxFile(fasta) as f:
        for kmerCounts in workers.imap_unordered(find_kmers, [ (str(Seq(entry.sequence, generic_dna).reverse_complement()), k, None) for entry in f]):
            finalKmerCounts+=kmerCounts
elif strand=='both':
    with pysam.FastxFile(fasta) as f:
        for kmerCounts in workers.imap_unordered(find_kmers, [ (str(entry.sequence), k, None) for entry in f]):
            finalKmerCounts+=kmerCounts
    with pysam.FastxFile(fasta) as f:
        for kmerCounts in workers.imap_unordered(find_kmers, [ (str(Seq(entry.sequence, generic_dna).reverse_complement()), k, None) for entry in f]):
            finalKmerCounts+=kmerCounts
else:
    print("Wrong strand specification (use '+' or '-' or 'both')")

# print("kmer\tabundance")
for kmer,abundance in finalKmerCounts.most_common():
    print(f"{kmer},{abundance}")
