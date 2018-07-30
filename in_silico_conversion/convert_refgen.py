import argparse
import copy
import itertools as it
from Bio.Seq import Seq,MutableSeq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys
sys.stdout.flush() # To flush output to std out

argparser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description="In silico bisulfite convert genome based on reference methylome \n By Emma Dann",
    )
argparser.add_argument('refgen', type=str, help='Fasta file of reference genome')
argparser.add_argument('refMet', type=str, help='Cov2c file of reference methylome ')
argparser.add_argument('--outputPrefix', default='mm10', type=str,help='Format of output file: bedGraph or bigWig')
args = argparser.parse_args()

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
    cgn = 0
    totCg = len(open(metFile, 'r').readlines())
    with open(metFile, 'r') as f:
        for cgSite in f.readlines():
            chr,pos,strand,met = cgSite.rstrip().split('\t')
            cgn += 1
            progress = round((cgn/totCg)*100,4)
            # print(chroms['chr8'].seq[int(pos)-1])
            if strand==strandOI:
                # rec = [conv for conv in convertedSeqs.values() if conv.name==chr]
                if met=='1':
                    convertedSeqs[chr].seq[int(pos)-1] = 'C'
            if progress % 5 == 0:
                print('Processed ' + str(progress) + '%  of met file')
    return(convertedSeqs)

def convert_seqs(seqs, metFile, strand):
    convertedSeqs = copy.deepcopy(seqs)
    for chr,rec in convertedSeqs.items():
        print('--- Making complete conversion ---')
        rec.seq = complete_conversion(rec.seq)
        print('--- Making methylome based deconversion ---')
        metConvertedSeqs = methylome_deconversion(convertedSeqs, metFile, strandOI=strand)
        print('--- Saving ---')
    return(metConvertedSeqs)

# small_recs={}
# for chr,rec in crickChroms.items():
#     rec.seq=rec.seq[0:10]
#     small_recs[chr]=rec

def convert_refgen(refgen, metFile, outputPrefix):
    chroms = read_refgen(refgen)
    print('--- Building reverse strand ---')
    revChroms = make_crick_strand(chroms)
    with open(outprefix + '_' + refgen.split('/')[-1].rstrip('.fa') + '.BSconv.forward.fa','w') as out:
        SeqIO.write(convert_seqs(chroms,metFile, '+').values(), out, 'fasta')
    with open(outprefix + '_' + refgen.split('/')[-1].rstrip('.fa') + '.BSconv.reverse.fa','w') as out:
        SeqIO.write(convert_seqs(revChroms,metFile, '-').values(), out, 'fasta')

# refgen = '/hpc/hub_oudenaarden/edann/hexamers/kaester/test_refgen.fa'
# # metFile = '/hpc/hub_oudenaarden/edann/hexamers/kaester/met_extraction/merged_reference_CG.met.chr1'
# # refgen='../test_refgen.fa'
# metFile='/hpc/hub_oudenaarden/edann/hexamers/kaester/test_CG_met.txt'
refgen=args.refgen
metFile=args.refMet
outputPrefix=args.outputPrefix

convert_refgen(refgen, metFile, outputPrefix)
