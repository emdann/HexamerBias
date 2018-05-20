#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) bamfile"
    echo "2) reference genome (mm10, hg38, danRer10, WBcel235)"
    echo "3) untrimmed fasta"
    echo "4) type of mapping (bs_se or bs_pe)"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bamfile=$1
refgen=$2
fasta=$3
type=$4
sample=$(echo $bamfile | sed 's,/,\t,g' | awk '{gsub(/\..+/, ""); print $NF}')
genome=$(echo $refgen | sed 's,/,\t,g' | awk '{gsub(/\..+/, ""); print $NF}')
abundanceDir='/hpc/hub_oudenaarden/edann/hexamers/genome_kmers'
abfile=${abundance_dir}/${genome}.kmerAbundance.csv

bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction

echo "--- Step 1: get primed region --- "
python ${bin_dir}/getPrimedRegion.py -o ./ -t $type -s $bamfile $refgen
echo "--- Step 2: Make pt table --- "
python ${bin_dir}/bsPrimerTemplTab.py $fasta ${sample}.primedreg.fa $abfile
echo "--- Step 3: predict delta G --- "
python ${bin_dir}/ptModel.py ${sample}.ptCounts.qualFilt.parallel.csv $abfile bs

rm ${sample}.primedreg.fa
rm ${sample}.primedreg.bed
gzip ${sample}.ptCounts.qualFilt.parallel.csv
gzip ${sample}_ptDg_qual.csv
