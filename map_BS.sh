#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) fastq.gz file"
    echo "2) reference genome folder"
    echo "3) Type of data (BS or noBS)"
    echo "3) output directory"
    exit
fi

file=$1
refgen=$2
type=$3
outdir=$4

sample=$(echo $file | awk '{gsub(/.fastq.gz/, ""); print}')

path_2_bismark=/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3
path_2_samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1
path_2_bowtie=/hpc/hub_oudenaarden/bin/software/bowtie2-2.2.9
path_2_cutadapt=/hpc/hub_oudenaarden/edann/venv2/bin/cutadapt
path_2_trimgalore=/hpc/hub_oudenaarden/edann/bin/TrimGalore-0.4.3
path_2_bwa=/hpc/hub_oudenaarden/bin/software/bwa-0.7.10

echo "---- Trimming! ----"
if [ -e ${sample}_trimmed.fq.gz ]
then
  echo "${sample}_trimmed.fq.gz found"
else
  echo "sample=$(echo $file | awk '{gsub(/.fastq.gz/, ""); print}'); ${path_2_trimgalore}/trim_galore --clip_R1 9 --three_prime_clip_R1 3 --path_to_cutadapt ${path_2_cutadapt} -o $outdir $file" | \
      qsubl -N trim_${sample}

echo "---- Mapping! ----"
if [["$type" == "BS"]]
then
  echo "${path_2_bismark}/bismark --samtools_path ${path_2_samtools} --path_to_bowtie  $refgen ${sample}_trimmed.fq.gz" | \
      qsub -cwd -N map_${sample} -pe threaded 10 -l h_rt=24:00:00 -l h_vmem=50G -l h_cpu=1:00:00 -hold_jid trim_${sample}
else
  echo "${path_2_bwa}/bwa mem -t 10 $refgen ${sample}_trimmed.fq.gz | ${path_2_samtools}/samtools view -bS - > ${sample}.bam" | \
      qsub -cwd -N map_${sample} -pe threaded 10 -l h_rt=24:00:00 -l h_vmem=50G -l h_cpu=1:00:00 -hold_jid trim_${sample}

echo "---- Deduplicating! ----"
echo "${path_2_bismark}/deduplicate_bismark --samtools_path ${path_2_samtools} -s --bam ${sample}_trimmed_bismark_bt2.bam" | //
    qsubl -N dedup_${sample} -hold_jid map_${sample}
