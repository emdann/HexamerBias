#!/bin/bash

# retrieve the random hexamers from BAM file identifiers


input=$1
path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1
input_bam=${input}_trim1_R1_bismark_bt2_pe.deduplicated.bam

## Sort and index bam file
$path2samtools/samtools sort --output-fmt BAM -T tmp_${input} -o sorted_$input_bam $input_bam
rm tmp_${input}*

## get fastq file of bam
#$path2samtools/samtools fastq sorted_$input_bam
## extract IDs for two read files
$path2samtools/samtools fastq sorted_$input_bam | grep -E "^@.*2$"| tr "_" " " | sed 's=/2==' > IDs_2_${input}.txt
$path2samtools/samtools fastq sorted_$input_bam | grep -E "^@.*1$"| tr "_" " " | sed 's=/1==' > IDs_1_${input}.txt
## extract hexamers
gunzip ${input}_R1.fastq.gz
gunzip ${input}_R2.fastq.gz

awk 'NR==FNR{a[$1];next} $1 in a{getline; print substr($0,0,6);}' IDs_2_${input}.txt ${input}_R2.fastq > ${input}_hex_2.txt
awk 'NR==FNR{a[$1];next} $1 in a{getline; print substr($0,0,6);}' IDs_1_${input}.txt ${input}_R1.fastq > ${input}_hex_1.txt

#cat ${input}_R2.fastq | grep -A1 -f  IDs_2.txt | grep -E '^[A-Z]' | awk '{print substr($N,0,6)}'> ${input}_hex_2.txt
#cat ${input}_R1.fastq | grep -A1 -f  IDs_1.txt | grep -E '^[A-Z]' | awk '{print substr($N,0,6)}'> ${input}_hex_1.txt

rm IDs_*_${input}.txt
