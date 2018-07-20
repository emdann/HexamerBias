#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give:"
    echo "1) bismark bam file"
    echo "2) reference genome (mm10, hg38, danRer10, WBcel235)"
    echo "3) output directory"
    exit
fi

bam=$1
refgen=$2
outdir=$3

# sample=$(echo $file | awk '{gsub(/.fastq.gz/, ""); print}')

path_2_bismark=/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3
path_2_samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1
path_2_bowtie=/hpc/hub_oudenaarden/bin/software/bowtie2-2.2.9
path_2_cutadapt=/hpc/hub_oudenaarden/edann/venv2/bin/cutadapt
path_2_trimgalore=/hpc/hub_oudenaarden/edann/bin/TrimGalore-0.4.3
path_2_bwa=/hpc/hub_oudenaarden/bin/software/bwa-0.7.10

## DEFINE THE REFERENCE GENOME
if [[ "$refgen" == "hg38" ]]
  then
    refgen_dir=/hpc/hub_oudenaarden/cgeisenberger/genomes/hg38
  elif [[ "$refgen" == "WBcel235" ]]
    then
      refgen_dir=/hpc/hub_oudenaarden/edann/genomes/WBcel235
  else
    refgen_dir=/hpc/hub_oudenaarden/avo/BS/${refgen}
fi
echo "Refgen folder: $refgen_dir"

## Extract methylation
echo "${path_2_bismark}/bismark_methylation_extractor --samtools_path ${path_2_samtools} -s --gzip --report --multicore 8 --comprehensive --merge_non_CpG --bedGraph -o . --genome_folder $refgen_dir $bam" | \\
  qsub -cwd -N meth_extract_${bam} -l h_rt=01:00:00 -l h_vmem=20G -l h_cpu=1:00:00
