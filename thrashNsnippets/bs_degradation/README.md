### BS MEDIATED DEGRADATION
Do I see a strand bias in the most covered regions? Are the most abundant hexamers always on a certain strand?

1) Take regions with high cov (>5)
```
bedtools genomecov -ibam ../crypts_bs/VAN1667/sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam -g ../genomes/mm10/mm10.genome -bg | awk '$4 >=5' > sorted_L1_trim1_R1_bismark_bt2_pe.highcov.bed
```
Putting in a plus strand information for all and merging the overlapping regions together:
```
cat sorted_L1_trim1_R1_bismark_bt2_pe.highcov.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"."$2"."$3"\t"$4"\t+"}'
```

Make sliding windows of a 100 bp size and select only windows with the set size
```
bedtools makewindows -b sorted_L1_trim1_R1_bismark_bt2_pe.highcov.strand.merge.bed -w 100 -s 5|
  awk '($3-$2)==100' |
  awk '{print $1"\t"$2"\t"$3"\t"$1"."$2"."$3"\t.\t+"}' > sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.bed
```

2) Check coverage for + and - strand
Plus strand (bedtools multicov -s searches for reads in the same strand):
```
bedtools multicov -bams ../crypts_bs/VAN1667/sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam -bed sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.bed -s > sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.plus.bed
```
Minus strand (bedtools multicov -s searches for reads in the opposite strand):
```
bedtools multicov -bams ../crypts_bs/VAN1667/sorted_L1_trim1_R1_bismark_bt2_pe.deduplicated.bam -bed sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.bed -S > sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.minus.bed
```

3) Take plus and minus information for every region and compute total coverage and strand ratio. Add nucleotide content information.
```
  bedtools intersect -a sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.plus.bed -b sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.minus.bed -wb |
  awk 'BEGIN {FS="\t"}
     {if($14!=0){ratio=$7/$14}
     else{ratio=1000}}
     {print $1,$2,$3,$4, $7+$14,$6,log(ratio)}' |
  sed 's/ /\t/g' |
  bedtools nuc -fi ../genomes/mm10/mm10.fa -bed stdin > sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.ratio.bed
```
