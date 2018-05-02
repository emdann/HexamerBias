#!/bin/bash

if [ $# < 1 ]
then
    echo "Please, give:"
    echo "1) Fasta file of reference genome"
    echo "Optional: "
    echo "2) Number of threads to use (default = 8)"
    exit
fi

genome=$1
threads=${2:-8} 

organism=$(echo $refgen | sed 's,.*/,,g' | sed 's/.fa//')
bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias/utils

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
python ${bin_dir}/kmerCounter.py -t $threads | grep -v N > ${organism}.kmerAbundance.csv
