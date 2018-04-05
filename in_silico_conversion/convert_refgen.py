import pandas as pd
import itertools as it
from Bio.Seq import Seq,MutableSeq,Records
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
    complRecs = {}
    for chr,rec in seqs.items():
        rec.seq = rec.seq.complement()
        complRecs[chr] = rec
    return(complRecs)

def complete_conversion(seq):
    '''
    Convert all Cs in the refgen dictionary to Ts
    '''
    newseq = Seq(''.join(('T' if n=='C' else n for n in list(seq))), generic_dna)
    return(newseq)

def methylome_deconversion(seqs, metFile):
    '''
    Revert conversion based on reference genome: convert Ts to C if the position
    is methylated.
    '''
    convertedSeqs = []
    for rec in seqs:
        rec.seq = MutableSeq(str(rec.seq), generic_dna)
        convertedSeqs.append(rec)
    with open(metFile, 'r') as f:
        for cgSite in f.readlines():
            chr,pos,strand,met = cgSite.rstrip().split('\t')
            if met=='1':
                convertedSeqs[chr][int(pos)] = 'C'
    return(convertedSeqs)

def convert_seqs(seqs, metFile):
    convertedSeqs = {}
    for chr,rec in seqs.items():
        rec.seq = complete_conversion(rec.seq)
        convertedSeqs[chr] = rec
    metConvertedSeqs = methylome_deconversion(convertedSeqs, metFile)
    return(metConvertedSeqs)

# small_recs={}
# for chr,rec in crickChroms.items():
#     rec.seq=rec.seq[0:10]
#     small_recs[chr]=rec

def convert_refgen(refgen, metFile):
    chroms = read_refgen(refgen)
    revChroms = make_crick_strand(chroms)
    with open('mm10.crypts.BSconv.forward.fa','w') as out:
        SeqIO.write(convert_seqs(chroms,metFile), out, 'fasta')
    with open('mm10.crypts.BSconv.reverse.fa','w') as out:
        SeqIO.write(convert_seqs(revChroms,metFile), out, 'fasta')

refgen = '/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa'
metFile = 'test_reference_CG.txt'
convert_refgen(refgen, metFile)
