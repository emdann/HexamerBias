##### IS METHYLATION STATUS DRIVING THE COVERAGE BIAS?
A number of reasons to point in this direction:
1) Melting temperature of the protocol is enough to open up unmethylated regions but not methylated ones
2) Unmethylated regions are more cut and so more covered

### Methylation levels in Reik compared to normal BS seq
Is the W.Reik protocol generally underestimating methylation?

- Take CpG methylation of kaester file and reik file (cov > 5):
```
zcat ERR454965_1_val_1_bismark_bt2.deduplicated.bismark.cov.gz | awk '$5+$6>=5' > ERR454965_1_val_1_bismark_bt2.deduplicated.bismark.thresh5.cov
zcat crypts_bs/VAN1667/L1_trim1_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '$5+$6>=5' >
crypts_bs/VAN1667/L1_trim1_R1_bismark_bt2_pe.deduplicated.bismark.thresh5.cov
```
