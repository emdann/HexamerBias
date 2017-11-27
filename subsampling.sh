#!/bin/bash 

smp=$1
p=$2

samtools view -b -s $p $smp > sub_${p}_${smp}
/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3/bismark_methylation_extractor --samtools_path /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 --comprehensive -s --genome_folder /hpc/hub_oudenaarden/avo/BS/mm10 sub_${p}_${smp}
rm sub_${p}_${smp}

