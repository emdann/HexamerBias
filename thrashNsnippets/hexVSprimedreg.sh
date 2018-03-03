#!/bin/sh

#  hexVSprimedreg.sh
#  
#
#  Created by Emma on 12/01/2018.
#
source /hpc/hub_oudenaarden/venv2/bin/activate

echo "source /hpc/hub_oudenaarden/edann/venv2/bin/activate; python ../bin/coverage_bias/getPrimedRegion.py ../crypts_bs/VAN1667/sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam > L1_primed_reg.qname.read.bed"

cat L1_primed_reg.qname.read.bed | awk '\$5==2' | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -fi ../genomes/mm10/mm10.fa -bed stdin -fo L1_R2_primed_seq.original.fa
cat L1_primed_reg.qname.read.bed | awk '\$5==1' | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools getfasta -name -fi ../genomes/mm10/mm10.fa -bed stdin -fo L1_R1_primed_seq.original.fa

python /hpc/hub_oudenaarden/bin/coverage_bias/hexVSprimed.py /hpc/hub_oudenaarden/edann/crypts_bs/VAN1667/L1_R2.fastq.gz /hpc/hub_oudenaarden/edann/hexamers/L1_R2_primed_seq.original.fa


