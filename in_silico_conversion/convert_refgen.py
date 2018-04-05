import pandas as pd
import copy
import itertools as it
from Bio.Seq import Seq,MutableSeq,Records
from Bio import SeqIO
from Bio.Alphabet import generic_dna

# Open refgen and convert everything
def read_refgen(refgen):
    seqs = SeqIO.to_dict(SeqIO.parse(refgen,'fasta'))
    upperSeqs = {}
    for chr,rec in seqs.items():
        rec.seq = rec.seq.upper()
        upperSeqs[chr] = rec
# # with ps.FastxFile(refgen) as fa:
    #  	# for entry in fa:
    #  		seqs[entry.name]= entry.seq.upper()
    return(upperSeqs)

def make_crick_strand(seqs):
    '''
    Makes list of records of complementary sequences (direction 3'->5'!!!)
    '''
    complRecs = copy.deepcopy(seqs)
    for chr,rec in complRecs.items():
        rec.seq = rec.seq.complement()
        # complRecs[chr] = rec
    return(complRecs)

def complete_conversion(seq):
    '''
    Convert all Cs in the refgen dictionary to Ts
    '''
    newseq = Seq(''.join(('T' if n=='C' else n for n in list(seq))), generic_dna)
    return(newseq)

def methylome_deconversion(seqs, metFile, strandOI='+'):
    '''
    Revert conversion based on reference genome: convert Ts to C if the position
    is methylated.
    '''
    mutableSeqs = copy.deepcopy(seqs)
    for chr,rec in mutableSeqs.items():
        rec.seq = MutableSeq(str(rec.seq), generic_dna)
        # mutableSeqs[chr]=rec
    convertedSeqs = mutableSeqs
    with open(metFile, 'r') as f:
        for cgSite in f.readlines():
            chr,pos,strand,met = cgSite.rstrip().split('\t')
            if strand==strandOI:
                # rec = [conv for conv in convertedSeqs.values() if conv.name==chr]
                if met=='1':
                    convertedSeqs[chr].seq[int(pos)] = 'C'
    return(convertedSeqs)

def convert_seqs(seqs, metFile, strand):
    convertedSeqs = copy.deepcopy(seqs)
    for chr,rec in convertedSeqs.items():
        rec.seq = complete_conversion(rec.seq)
        metConvertedSeqs = methylome_deconversion(convertedSeqs, metFile, strandOI=strand)
    return(metConvertedSeqs)

# small_recs={}
# for chr,rec in crickChroms.items():
#     rec.seq=rec.seq[0:10]
#     small_recs[chr]=rec

def convert_refgen(refgen, metFile):
    chroms = read_refgen(refgen)
    revChroms = make_crick_strand(chroms)
    with open('mm10.crypts.BSconv.forward.fa','w') as out:
        SeqIO.write(convert_seqs(chroms,metFile, '+').values(), out, 'fasta')
    with open('mm10.crypts.BSconv.reverse.fa','w') as out:
        SeqIO.write(convert_seqs(revChroms,metFile, '-').values(), out, 'fasta')

# refgen = '/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'
# metFile = 'test_reference_CG.txt'
refgen='../test_refgen.fa'
metFile='../test_CG_met.txt'

convert_refgen(refgen, metFile)
