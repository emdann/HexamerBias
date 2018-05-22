#!/bin/bash

bam=$1
sample=$(echo $bam | awk '{gsub(/.bam/, ""); print}')
path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1
sort_order=$(${path2samtools}/samtools view -H $bam | head -1 | cut -f 3 )

if [[ $sort_order =~ 'unsorted' ]] ;
then
  echo "Sorting file into ${sample}.srt.bam"
  ${path2samtools}/samtools sort -@ 5  ${bam} -o ${sample}.srt.bam
  bam=${sample}.srt.bam
else
  echo "Bamfile is sorted"
fi

echo "Indexing file ${sample}.srt.bam"
${path2samtools}/samtools index -@ 5 $bam
