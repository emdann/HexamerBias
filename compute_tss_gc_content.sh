#!/bin/bash
org=$1
flank=$2

slop=$flank/2

cat hgTables.refgen.${org}.txt | 
	awk '{print $3,$5,$5}' | 
	tr ' ' '\t' | 
	tail -n+2 | 
	grep -v '_'| 
	/hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools slop -b $slop -i stdin -g ../genomes/${org}/${org}.genome | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools nuc -fi ../genomes/${org}/${org}.fa -bed stdin | cut -f 1,2,3,5 > GCcont_txEnd_${flank}slop_hgTables.refgen.${org}.txt
