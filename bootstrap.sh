#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Please, give:"
    echo "1) bamfile"
    echo "2) reference genome"
    echo "3) untrimmed fasta"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
bamfile=$1
refgen=$2
fasta=$3
abfile=/hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz
sample=$(echo $bamfile | awk '{gsub(/.bam/, ""); print}')

if [ $refgen = "mouse" ]
then
    refgen=/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa
fi

bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias

echo "--- Step 1: get primed region --- "
python ${bin_dir}/getPrimedRegion.py -o ./ $bamfile $refgen bs_se
echo "--- Step 2: Make pt table --- "
python ${bin_dir}/bsPrimerTemplTab.py $fasta ${sample}.primedreg.fa $abfile
echo "--- Step 3: predict delta G --- "
python ${bin_dir}/ptModel.py ${sample}.ptCounts.qualFilt.parallel.csv $abfile bs

rm ${sample}.primedreg.fa
rm ${sample}.primedreg.bed
gzip ${sample}.ptCounts.qualFilt.parallel.csv
gzip ${sample}_ptDg_qual.csv
