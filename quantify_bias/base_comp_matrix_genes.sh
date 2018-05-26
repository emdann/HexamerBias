#!/bin/bash

gen_fasta=/hpc/hub_oudenaarden/gene_models/cel_gene_models/WBcel235.fa
bed=/hpc/hub_oudenaarden/edann/hexamers/annotations_bed/RefSeq_genes_WBcel235_noChr.bed
gen_file=/hpc/hub_oudenaarden/edann/genomes/WBcel235.genome
bins=$1
slop=$2

cat $bed  | cut -f 1,2,3,4,5,6 | grep '+' | bedtools makewindows -b stdin -n $bins -i winnum > genes_plus.bed
cat $bed  | cut -f 1,2,3,4,5,6 | grep '-' | bedtools makewindows -reverse -b stdin -n $bins -i winnum > genes_minus.bed

cat genes_*us.bed | \
  bedtools nuc -fi $gen_fasta -bed stdin | \
  sort -k4n | \
  bedtools groupby -i stdin -g 4 -c 5,6 -o mean,mean > nuc_refseq_genes_WBcel235.txt

rm genes_*us.bed
