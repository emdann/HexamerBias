#!/bin/bash
genes_bed=$1
refgen_fa=$2
genome=$3



path2bedtools=/hpc/hub_oudenaarden/edann/bin/bedtools2/bin

cat $genes_bed | \
	cut -f 1,2,3 | \
	awk '$3=$2{print}' | \
	grep -v '_' |
	tr ' ' '\t' > TxStart_${genes_bed}

cat TxStart_${genes_bed} | \
	${path2bedtools}/bedtools slop -b 50 -i stdin -g $genome | \
	${path2bedtools}/bedtools nuc -fi $refgen_fa -bed stdin | \
	cut -f 1,2,3,5 > TxStart_${genes_bed}.GCcont.bed

	cat TxStart_${genes_bed} | \
		${path2bedtools}/bedtools slop -b $slop -i stdin -g $genome | \
		${path2bedtools}/bedtools shuffle -chrom -i stdin -g $genome | \
		${path2bedtools}/bedtools nuc -fi $refgen_fa -bed stdin | \
		cut -f 1,2,3,5 > TxStart_shuffle_${genes_bed}.GCcont.bed
