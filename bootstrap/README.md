### BOOTSTRAP ON POOLED SAMPLES

1) Pool all samples from the same experiment in one bam file
```
samtools cat crypts_bs/VAN1667/sorted_L*.bam > VAN1667.bam
```
2) Make subsamples of different percentages of reads of the pooled bam (run on multiple cores)
```
for f in $(seq 15 5 95);
  do
  echo "/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools view -s 0.$f  -o VAN1667.$f.bam -@ 8 VAN1667.bam" |
  qsubl -N subsmp.${f} -pe threaded 9;
  done
```

3) Run DeltaF prediction on all subsamples
```
for file in VAN1667_[123456789]*.bam;
  do
  echo "/hpc/hub_oudenaarden/edann/bin/coverage_bias/run_deltaF_prediction.sh $file mouse VAN1667.fastq.gz" |
  qsubl -N predDg_${file}; done
```

4) Make table of predicted coverage for each predicted DeltaF table
```

```