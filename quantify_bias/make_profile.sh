#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) bamfile"
    echo "2) regions of interest (bed or gtf file)"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bamfile=$1
bed=$2
sample=$(echo $bamfile | awk '{gsub(/.bam/, ""); print}')
regions=$(echo $bed | awk '{gsub(/.bed/, ""); print}')
## Make coverage BW
echo "--- Computing coverage ---"
bamCoverage -b $bamfile -o ${sample}.bw

## Make coverage matrix on defined region
echo "--- Computing matrix ---"
computeMatrix scale-regions \
  -R  $bed \
  -S ${sample}.bw  \
  -b 3000 -a 3000 \
  --regionBodyLength 5000 --skipZeros \
  -o ${sample}.mat.gz \
  --outFileNameMatrix ${sample}.mat.tab

## Plot profile
echo "--- Plotting profile ---"
plotProfile -m ${sample}.mat.gz \
              -out ${sample}_${regions}_coverage.png  \
              --numPlotsPerRow 1 \
              --yAxisLabel "coverage" \
              --samplesLabel $sample
