#!/bin/bash

input=$1

/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools view ${input}_trim1_R1*deduplicated.bam | awk '$2==163 || $2==83 ' | cut -f 1 | uniq > OB_${input}.test
/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools view ${input}_trim1_R1*deduplicated.bam | awk '$2==99 || $2==147' | cut -f 1 | uniq > OT_${input}.test

cat hexN10bps_${input}.txt | tr '/' '\t' | tr -d '>' > hexN10bps_${input}.test
awk 'NR==FNR{a[$1];next} $1 in a{if ($2==1) {print $N"\tOT"} else {print $N"\tCTOT"}}' OT_${input}.test hexN10bps_${input}.test > hexN10bps_${input}_OT.txt 
awk 'NR==FNR{a[$1];next} $1 in a{if ($2==1) {print $N"\tOB"} else {print $N"\tCTOB"}}' OB_${input}.test hexN10bps_${input}.test > hexN10bps_${input}_OB.txt

cat hexN10bps_${input}_OT.txt hexN10bps_${input}_OB.txt > hexN10bps_str_${input}.txt

#rm hexN10bps_${input}_O*.txt
rm *${input}.test


