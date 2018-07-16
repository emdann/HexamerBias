#!/bin/bash

## Compute total genomic coverage ##

bamfile=$1
path2bedtools=/hpc/hub_oudenaarden/edann/bin/bedtools2/bin

${path2bedtools}/bedtools genomecov -ibam ${bamfile} | awk '\$2==0' | grep genome | cut -f 5
