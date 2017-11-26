#!/bin/bash

zcat refGen.hg38.txt.gz | awk '{print $3,$5,$5}' | tr ' ' '\t' | tail -n+2 | grep -v '_'| /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools slop -b 500 -i stdin -g ../genomes/hg38.genome | /hpc/hub_oudenaarden/edann/bin/bedtools2/bin/bedtools nuc -fi ../genomes/hg38.fa -bed stdin | cut -f 1,2,3,5 > GCcont_1000slop_hgTables.refgen.hg38.txt