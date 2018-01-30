#!/bin/sh

#  findSNPsPerPosition.sh
#
#
#  Created by Emma on 10/01/2018.
#
if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) input file (.bam,.sam)"
    echo "2) ref genome (path)"
    exit
fi

path2samtools=/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1

bam=$1
refgen=$2
name=$(echo $bam | sed 's/..am$//')

${path2samtools}/samtools mpileup -O -I -f $refgen $bam > $name.mpileup

cat $name.mpileup |
awk '$5 ~ /[ACGTNacgtn]/ {print}' | # Select lines with mismatches
awk '!($3 ~ /[cC]/ && $5=="t" || $5=="T")' | # Rule out lines with C to T conversion
awk '!($3 ~ /[gG]/ && $5=="a" || $5=="A")' | # Rule out lines with G to A conversion
# awk '{if($7 ~ /,/){print }}' # Split lines with multiple reads
awk '{if($5 ~ /[\^\+\-]./){gsub(/[\^\+\-]./, "")}}{print}' | # Remove extra chars for start and 1 bp indel
sed 's/\$//'| # Remove extra chars for end
awk '
    {
    n = split($7, pos, ",")
    m = split($5, baseRead, "")
    p = split($6, qual, "")

    for (i=1; i<=n; i++)
        print $1, $2, $3, $4, baseRead[i], qual[i], pos[i]
    }
    ' |
tr ' ' '\t' |
awk '$5 ~ /[ACGTNacgtn]/ {print}' |
awk '!($3 ~ /[cC]/ && $5=="t" || $5=="T")' |
awk '!($3 ~ /[gG]/ && $5=="a" || $5=="A")' > $name.pileup
