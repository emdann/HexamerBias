#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give:"
    echo "1) bamfile"
    echo "2) reference genome fasta file (mm10, hg38, danRer10, WBcel235)"
    echo "3) untrimmed fasta"
    echo "4) type of mapping (bs_se or bs_pe or no_bs)"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bamfile=$1
refgen=$2
fasta=$3
type=$4
sample=$(echo $bamfile | sed 's,/,\t,g' | awk '{gsub(/\..+/, ""); print $NF}')
genome=$(echo $refgen | sed 's,/,\t,g' | awk '{gsub(/\..+/, ""); print $NF}')
abundance_dir='/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers'
abfile=${abundance_dir}/${genome}.kmerAbundance.csv

bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction

if [ -e ${sample}.primedreg.fa ]
then
    echo "Fasta of template regions found!"
else
  echo "--- Step 1: get primed region --- "
  python ${bin_dir}/getPrimedRegion.py -o ./ -t $type -s $bamfile $refgen
fi

echo "--- Step 2: Make pt table --- "
python ${bin_dir}/bsPrimerTemplTab.py -t $type $fasta ${sample}.primedreg.fa $abfile

# echo "--- Step 3: predict delta G --- "
# python ${bin_dir}/ptModel.py ${sample}.ptCounts.qualFilt.csv $abfile

rm ${sample}.primedreg.fa
rm ${sample}.primedreg.bed
# gzip ${sample}.ptCounts.qualFilt.parallel.csv
# gzip ${sample}_ptDg_qual.csv
