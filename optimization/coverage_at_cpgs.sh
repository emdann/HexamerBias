#!/bin/bash

bw_file=$1
sample=$(echo $bw_file | sed 's/.bw//')
cpgs_file=/hpc/hub_oudenaarden/edann/hexamers/annotations_bed/cpg_sites_mm10.bed

bigWigToWig ${sample}.bw ${sample}.wig
convert2bed --input=wig < ${sample}.wig > ${sample}.bed

cat ${sample}.bed | \
  bedtools intersect -a stdin -b $cpgs_file > ${sample}.cpg_coverage.bed

n_cpgs=$(wc -l ${sample}.cpg_coverage.bed | cut -f 1 -d ' ')
bedtools sample -n $n_cpgs -i ${sample}.bed > ${sample}.random.bed
