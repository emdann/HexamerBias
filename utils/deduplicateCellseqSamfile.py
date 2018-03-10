import pysam as ps
import collections
import argparse

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Deduplicate PCR duplicated from CellSeq bamfile \n By Emma Dann")
argparser.add_argument('bamfile', type=str, help='Bam to deduplicate')
args = argparser.parse_args()

def findReadsToKeep(bamfile):
    '''
    Find list of reads to keep in deduplicated files.
    If more than one read with same position, umi, cell, primer, keep the first one
    '''
    bamDic = {}
    with ps.AlignmentFile(bamfile,"rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.flag==0:
                umi = r.qname.split(':')[-3]
                transcript = bam.getrname(r.reference_id)
                pos = r.pos
                primer = r.seq[0:6]
                if (umi, transcript, pos, primer) not in bamDic.keys():
                    bamDic[(umi, transcript, pos, primer)]=[]
                bamDic[(umi, transcript, pos, primer)].append(r.qname)
    # If more than one read with same position, umi, cell, primer, keep the first one
    toKeep = [names[0] for names in bamDic.values()]
    return(toKeep)

def saveDeduplicatedBam(bamfile, readsToKeep):
    '''
    Saves deduplicated version of bam file removing reads that have the same
    starting position, umi, cell, primer sequence.
    '''
    with ps.AlignmentFile(bamfile,"rb") as bam:
        new_file = ps.Samfile(bamfile.split('.')[0]+'.deduplicated.bam', mode='wb', template=bam)
        for r in bam.fetch(until_eof=True):
            if r.qname in toKeep:
                 new_file.write(r)
    return("Deduplicated bamfile saved.")

readsToKeep = findReadsToKeep(bamfile)
saveDeduplicatedBam(bamfile, readsToKeep)
