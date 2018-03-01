import pysam as ps
import collections

bamfile='/hpc/hub_oudenaarden/edann/hexamers/rnaseq/mouse/testing/SvdB11d1-MitoTrackerThird-Satellites-Adult_sorted.bam'

bamDic = {}
with ps.AlignmentFile(bamfile,"rb") as bam:
    for r in bam.fetch(until_eof=True):
        if r.flag==0:
            umi = r.qname.split(':')[-3]
            transcript = bam.getrname(r.reference_id)
            pos = r.pos
            if (umi, transcript, pos) not in bamDic.keys():
                bamDic[(umi, transcript, pos)]=[]
            bamDic[(umi, transcript, pos)].append(r.qname)

dupItems=[(loc,names) for loc,names in bamDic.items() if len(names)>1]

for loc,names in dupItems:
    umi,transcript,pos = loc
with ps.AlignmentFile(bamfile,"rb") as bam:
    for r in bam.fetch(transcript, start=pos):
        if r.flag==0 and sum(r.query_qualities[0:6]) >=196:
            print(r.seq[0:6], sum(r.query_qualities[0:6]))
                umi1 = r.qname.split(':')[-3]
                transcript1 = bam.getrname(r.reference_id)
                pos1 = r.pos

                if (umi1,transcript1,pos1) == (umi,transcript,pos):
                    print(r.seq[0:6])




def write_deduplicated_bam(bamfile):
    header = bamfile.header.copy()
    out = pysam.Samfile(outfile, 'wb', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)
