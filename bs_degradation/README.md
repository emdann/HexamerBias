### BS MEDIATED DEGRADATION
Do I see a strand bias in the most covered regions? Are the most abundant hexamers always on a certain strand?

1) Take regions with high cov (>5)
```
zcat VAN1667_depth.bed.gz | awk '$4>=5' > VAN1667_depth.highcov.bed
```
Putting in a 'pseudo strand' information:
```
cat VAN1667_depth.highcov.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"."$2"."$3"\t"$4"\t+"}' > VAN1667_depth.highcov.strand.bed
```

2) Check coverage for + and - strand
Same strand:
```
bedtools multicov -bams ../crypts_bs/VAN1667/sorted_L?_trim1_R1_bismark_bt2_pe.deduplicated.bam -bed VAN1667_depth.highcov.strand.bed -s > VAN1667.highcov.plusstrand.bed"
```
Opposite strand:
```
bedtools multicov -bams ../crypts_bs/VAN1667/sorted_L?_trim1_R1_bismark_bt2_pe.deduplicated.bam -bed VAN1667_depth.highcov.strand.bed -S > VAN1667.highcov.minusstrand.bed"
```
