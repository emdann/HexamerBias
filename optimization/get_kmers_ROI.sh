#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) bedfile of regions of interest (move to that folder!!!)"
    echo "2) folder to genome fasta files"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bed=$1
genomedir=$2
sample=$(echo $bed | awk '{gsub(/.bed/, ""); print}')
genome=$(echo $genomedir | awk -F"/" '{print $NF}')

bindir=/hpc/hub_oudenaarden/edann/bin/coverage_bias
path2bedtools=/hpc/hub_oudenaarden/edann/bin/bedtools2/bin

## Kmer count ROI

${path2bedtools}/bedtools getfasta -fi ${genomedir}/${genome}.fa -bed $bed -fo ${sample}.fa
${bindir}/utils/kmerCounter.py -s both ${sample}.fa > ${sample}.kmersTot.csv

## Kmer count random regions
${path2bedtools}/bedtools shuffle -seed 42 -chrom  -i $bed -g ${genomedir}/${genome}.genome | \
${path2bedtools}/bedtools getfasta -fi ${genomedir}/${genome}.fa -bed stdin -fo ${sample}.randomize.fa
${bindir}/utils/kmerCounter.py -s both ${sample}.randomize.fa > ${sample}.randomize.kmersTot.csv

## FoldChange
${bindir}/optimization/compute_kmers_FC.py ${sample}.kmersTot.csv ${sample}.randomize.kmersTot.csv
grep -v 'N' ${sample}.kmersFC.csv > ${sample}.noN.kmersFC.csv
