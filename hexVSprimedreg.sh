#!/bin/sh

#  hexVSprimedreg.sh
#  
#
#  Created by Emma on 12/01/2018.
#

echo "source /hpc/hub_oudenaarden/edann/venv2/bin/activate; python ../bin/coverage_bias/getPrimedRegion.py ../crypts_bs/VAN1667/sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam > L1_primed_reg.qname.read.bed"

echo "cat L1_primed_reg.qname.read.bed | awk '\$5==2' | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -fi ../genomes/mm10/mm10.fa -bed stdin -fo L1_R2_primed_seq.original.fa" | qsub -cwd -N primedseq -l h_rt=10:00:00

echo "cat L1_primed_reg.qname.read.bed | awk '\$5==1' | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -fi ../genomes/mm10/mm10.fa -bed stdin -fo L1_R1_primed_seq.original.fa" | qsub -cwd -N primedseq -l h_rt=10:00:00
