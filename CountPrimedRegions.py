import itertools as it
import pandas as pd
import collections
import numpy as np
import pysam as ps
import sys
import os
import argparse
import multiprocessing
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

fasta='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/gk2a-2_primed_seq.fa'

templDic={}
with ps.FastxFile(fasta) as templ:
	for entry in templ:
		seq, name = entry.sequence.upper(), entry.name
		templDic[name]=seq
