#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) reference genome (path or mouse/zebrafish/zebrafishGFP)"
    echo "2) root for fastq files"
    exit
fi

# input sorted file 
input=$1
#refGen=$1

if [ $refGen = 'mouse' ]
then
    refGen=/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa
elif [ $refGen = 'human' ]
then
    refGen=/hpc/hub_oudenaarden/edann/genomes/hg19/hg19.fa
elif [ $refGen = 'zebrafishGFP' ]
then
    refGen=/hpc/hub_oudenaarden/aalemany/refSeq/BS-zf-lintrace_GFP
fi
path2bedtools=/hpc/hub_oudenaarden/edann/bin/bedtools2/bin


## Making fasta file of 10 bps upstream of mapped read in reference genome
#../../bin/bedtools2/bin/bedtools bamtobed -i sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam | cut -f 1,2,3 | ../../bin/bedtools2/bin/bedtools flank -l 10 -r 0 -g /hpc/hub_oudenaarden/edann/hg19/* -i stdin
#$path2bedtools/bedtools bamtobed -i sorted_${input}_trim1_R1_bismark_bt2_pe.deduplicated.bam | cut -f 1,2,3,4 | awk '{print $1"\t"$2-10"\t"$2"\t"$4}' > ${input}_10bps_upstream.bed
#$path2bedtools/bedtools getfasta -name -fi $refGen -bed ${input}_10bps_upstream.bed -fo  sorted_${input}_trim1_R1_bismark_bt2_pe.deduplicated_10bps_upstream.fasta

#rm ${input}_10bps_upstream.bed 

## Computing hexamers fasta file
cat sorted_${input}_trim1_R1_bismark_bt2_pe.deduplicated_10bps_upstream.fasta | grep -e '^>' | tr '_' ' ' | grep -e '1$'| sed 's!/1!!'| sed 's/>/@/' > IDs_1_${input}.txt
awk 'NR==FNR{a[$1];next} $1 in a{getline seq; print  $N"\t"substr(seq,0,6);}' IDs_1_${input}.txt ${input}_R1.fastq | tr ' ' '_' | tr '@' '>'| awk '{print $1"/1\t"$2}' > hex_1_${input}.txt
cat sorted_${input}_trim1_R1_bismark_bt2_pe.deduplicated_10bps_upstream.fasta | grep -e '^>' | tr '_' ' ' | grep -e '2$'| sed 's!/2!!'| sed 's/>/@/' > IDs_2_${input}.txt
awk 'NR==FNR{a[$1];next} $1 in a{getline seq; print  $N"\t"substr(seq,0,6);}' IDs_2_${input}.txt ${input}_R2.fastq | tr ' ' '_' | tr '@' '>'| sed 's/_2:/_1:/'| awk '{print $1"/2\t"$2}' > hex_2_${input}.txt
cat hex_1_${input}.txt hex_2_${input}.txt > hex_${input}.txt

cat sorted_${input}_trim1_R1_bismark_bt2_pe.deduplicated_10bps_upstream.fasta | awk '{getline line; print $N"\t"line}' > tmp_10bps_${input}.txt

awk 'FNR==NR{a[$1]=$2;next}{print $0,a[$1]}' hex_${input}.txt tmp_10bps_${input}.txt | tr [:lower:] [:upper:] > hexN10bps_${input}.txt
