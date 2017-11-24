#!/bin/bash

cat hgTables.refgen.mm10.txt | awk '{print $3,$5,$5}' | tr ' ' '\t' | tail -n+2 | grep -v '_'| /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools slop -b 500 -i stdin -g ../genomes/mm10/mm10.genome | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools nuc -fi ../genomes/mm10/mm10.fa -bed stdin | cut -f 1,2,3,5 > GCcont_1000slop_hgTables.refgen.dr10.txt