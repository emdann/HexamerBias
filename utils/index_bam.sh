#!/bin/bash

bam=$1
sample=$(echo $bamfile | awk '{gsub(/.bam/, ""); print}')
path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1

${path2samtools}/samtools sort -@ 5  ${bam} -o ${sample}.srt.bam
${path2samtools}/samtools index -@ 5 ${sample}.srt.bam
