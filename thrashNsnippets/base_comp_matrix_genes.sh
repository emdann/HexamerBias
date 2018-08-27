#!/bin/bash

gen_fasta=/hpc/hub_oudenaarden/gene_models/zebrafish_gene_models/danRer10_clean.fa
bed=/hpc/hub_oudenaarden/edann/hexamers/annotations_bed/RefSeq_genes_danRer10.bed
gen_file=/hpc/hub_oudenaarden/edann/genomes/danRer10.genome
bins=$1
slop=$2

bedtools sample -n 10000 -i $bed > smp_bed.bed
bed=smp_bed.bed

cat $bed  | cut -f 1,2,3,4,5,6 | grep '+' | \
  bedtools makewindows -b stdin -n $bins -i winnum | \
  awk '{$4=$4+300; print}' > genes_plus.bed
cat $bed  | cut -f 1,2,3,4,5,6 | grep '-' | \
  bedtools makewindows -reverse -b stdin -n $bins -i winnum | \
  awk '{$4=$4+300; print}' > genes_minus.bed

# Left flanks
cat $bed  | cut -f 1,2,3,4,5,6 | \
  bedtools flank -s -i stdin -g $gen_file -l 3000 -r 0| \
  bedtools makewindows -b stdin -n 300 -i winnum > genes_left_flank.bed
# Right flanks
cat $bed  | cut -f 1,2,3,4,5,6 | \
  bedtools flank -s -i stdin -g $gen_file -l 0 -r 3000| \
  awk '$1 !~ /_/' | \
  bedtools makewindows -b stdin -n 300 -i winnum | \
  awk '{$4=$4+800; print}' > genes_right_flank.bed

cat genes_*.bed | \
  tr ' ' '\t' | \
  bedtools nuc -fi $gen_fasta -bed stdin | \
  tail -n+2 | \
  sort -k4n | \
  bedtools groupby -i stdin -g 4 -c 5,6 -o mean,mean > nuc_refseq_genes_danRer10smp.txt

rm genes_*us.bed
