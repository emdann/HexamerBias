#!/bin/bash

bam=$1
path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/bin/

smp=$(echo $bam | sed 's/trim1_//' | sed 's/_bismark_bt2.deduplicated.bam//')
#echo $smp
gunzip ${smp}.fastq.gz

$path2samtools/samtools view $bam | awk '{if ($15=="XR:Z:CT"){print $1}else{next}}'| tr '_' ' ' | sed 's/^/@/' > original_${smp}.ids
awk 'NR==FNR{a[$1];next} $1 in a{getline line; print line}' original_${smp}.ids ${smp}.fastq > original_${smp}.fasta

$path2samtools/samtools view $bam | awk '{if ($15=="XR:Z:GA"){print $1}else{next}}'| tr '_' ' ' | sed 's/^/@/' > compl_${smp}.ids
awk 'NR==FNR{a[$1];next} $1 in a{getline line; print line}' compl_${smp}.ids ${smp}.fastq > compl_${smp}.fasta

rm *${smp}.ids
gzip ${smp}.fastq
