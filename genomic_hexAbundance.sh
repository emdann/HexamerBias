#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) reference genome (mouse)"
    echo "2) root to cov2c file"
    exit
fi

source /hpc/hub_oudenaarden/edann/venv2/bin/activate
genome=$1
smp=$2

if [ $genome = "mouse" ]
then
    genome=/hpc/hub_oudenaarden/edann/genomes/mm10
fi

chromosomes=$(ls /hpc/hub_oudenaarden/edann/genomes/mm10/chr*.fa.gz | sed 's/_[GUJ]/\t/' | cut -f 1 | grep .fa.gz | sed 's,.*/,,' | sed 's/.fa.gz//')
python /hpc/hub_oudenaarden/edann/bin/BS_binning/splitBychr.py ${smp}.gz ${smp}.

for chr in $chromosomes
	do
		zcat ${smp}.${chr} | awk '$4+$5!=0{print $N"\t"$4/($4+$5)}' > ${smp}.${chr}.frac
		echo "python /hpc/hub_oudenaarden/edann/bin/weighted_kmers.py ${genome}/${chr}.fa.gz ${smp}.${chr}.frac > ${smp}.${chr}.hex" #|
		# qsub -cwd -N hex_${chr}_${smp} -l h_vmem=50G -l h_rt=10:00:00 -pe threaded 3
	done
	


