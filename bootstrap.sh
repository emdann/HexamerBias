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

if [ $genome = "mouse" ]
then
    genome=/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa
fi

bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias

python ${bin_dir}/getPrimedRegion.py $bamfile $ref_gen bs_se -o ./
python ${bin_dir}/bsPrimerTemplTab.py $fasta ${sample}.primedreg.fa $abfile
python ${bin_dir}/ptModel.py ${sample}.ptCounts.qualFilt.parallel.csv $abfile bs
