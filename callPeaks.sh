#!/bin/bash
# Find coverage peaks in scBS-seq bam file

if [ $# -ne 3 ]
then
    echo "Please, give:"
    echo "1) input file (.bam,.sam)"
    echo "2) output directory"
    echo "3) ref genome (mm,hg)"
    exit
fi

path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1

bam=$1
outdir=$2
gen=$3

format=$(echo $bam | awk '{if(/.bam$/){ print "BAM"}else{print "SAM"}}')
name=$(echo $bam | tr '/' '\t' | awk '{print $NF}' | sed 's/..am$//')

macs2 callpeak -t $bam --outdir $outdir -n $name --format $format -g $gen --nomodel --keep-dup all
cat ${name}_peaks.narrowPeak | awk '$7>4.5' > ${name}_peaks.filt.narrowPeak
