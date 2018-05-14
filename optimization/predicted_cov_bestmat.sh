#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) sample (output of DE optimization)"
    echo "2) Test bed regions"
    exit
fi

sample=$1
testbed=$2
bindir=/hpc/hub_oudenaarden/edann/bin/coverage_bias
venv=/hpc/hub_oudenaarden/edann/venv2/bin/activate

kmerAbFile=/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/mm10.kmerAbundance.csv
genomefa=/hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa

source $venv
${bindir}/optimization/coverage_best_matrix.py $sample
echo "source $venv; python ${bindir}/artificial_coverage/strand_specific_artificial_coverage.py $kmerAbFile ${sample}.bestMat.coverage.csv $testbed $genomefa -t 10" | \
  qsub -cwd -N predictedcov_${sample} -pe threaded 10 -l h_rt=24:00:00 -l h_vmem=80G -l h_cpu=10:00:00
