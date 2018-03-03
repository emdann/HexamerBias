#!/bin/bash


input=$1
path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1
input_bam=${input}_bismark_bt2.deduplicated.bam

#$path2samtools/samtools sort --output-fmt BAM -T tmp_${input} -o sorted_$input_bam $input_bam


#$apth2samtools/samtools fastq sorted_$input_bam | grep -E "^@.*2$"| tr "_" " " | sed 's=/2==' > IDs_2_${input}.txt
$path2samtools/samtools fastq sorted_$input_bam | grep -E "^@"| tr "_" " " | sed 's=/1==' > IDs_${input}.txt

fq_input=$(echo $input | sed 's/trim1_//')
#echo $fq_input
#awk 'NR==FNR{a[$1];next} $1 in a{getline line;getline;getline qual; print substr(qual,0,6)"\t"substr(line,0,6)}' IDs_2_${input}.txt ${input}_R2.fastq  | awk '$1 ~ /^[A-Z]+$/ {print $2}' > ${input}_2.hq.hex
gunzip ${fq_input}.fastq.gz
awk 'NR==FNR{a[$1];next} $1 in a{getline line;getline;getline qual; print substr(qual,0,6)"\t"substr(line,0,6)}' IDs_${input}.txt ${fq_input}.fastq  | awk '$1 ~ /^[A-Z]+$/ {print $2}' > ${input}.hq.hex

gzip ${fq_input}.fastq
