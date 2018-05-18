#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Please, give:"
  echo "1) bamfile (indexed)"
  echo "2) regions of interest (bed or gtf file)"
  exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bamfile=$1
bed=$2
sample=$(echo $bamfile | awk '{gsub(/.bam/, ""); print}')
regions=$(echo $bed | awk '{gsub(/.bed/, ""); print}')

## Index bam file
if [ -e ${bamfile}.bai ]
then
    echo "Index file ${bamfile}.bai found"
else
    echo "--- Indexing bam file ---"
    /hpc/hub_oudenaarden/edann/bin/coverage_bias/utils/index_bam.sh $bamfile
    bamfile=$sample.srt.bam
    sample=$(echo $bamfile | awk '{gsub(/.bam/, ""); print}')
fi

## Make coverage BW
if [ -e ${sample}.bw ]
then
  echo "Coverage file ${sample}.bw found"
else
  echo "--- Computing coverage ---"
  bamCoverage -b $bamfile -o ${sample}.bw -p 6
fi

## Make coverage matrix on defined region
echo "--- Computing matrix ---"
computeMatrix scale-regions \
  -R  $bed \
  -S ${sample}.bw  \
  -b 3000 -a 3000 \
  --regionBodyLength 5000 --skipZeros \
  -p 6 \
  -o ${sample}.mat.gz \
  --outFileNameMatrix ${sample}.mat.tab

computeMatrix reference-point \
  -R  $bed \
  -S ${sample}.bw  \
  --referencePoint center \
  -b 3000 -a 3000 \
  --skipZeros \
  -p 6 \
  -o ${sample}.CTCF.mat.gz \
  --outFileNameMatrix ${sample}.CTCF.mat.tab

## Plot profile
echo "--- Plotting profile ---"
plotProfile -m ${sample}.mat.gz \
              -out ${sample}_${regions}_coverage.png  \
              --numPlotsPerRow 1 \
              --yAxisLabel "coverage" \
              --samplesLabel $sample
