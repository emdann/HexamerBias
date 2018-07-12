### BOOTSTRAP ON POOLED SAMPLES

1) Pool all samples from the same experiment in one bam file
```
samtools cat crypts_bs/VAN1667/sorted_L*.bam > VAN1667.bam
```
2) Make subsamples of different percentages of reads of the pooled bam (run on multiple cores)
```
for f in $(seq 15 5 95);
  do
  echo "/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools view -s 0.$f  -o ${sample}_$f.bam -@ 8 ${sample}.bam" |
  qsubl -N subsmp.${f} -pe threaded 9;
  done
```

3) Build pt table for all subsamples
```
for file in *.bam;
  do
  echo "/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/run_deltaF_prediction.sh $file $refgen $fasta $type" |
  qsubl -N predDg_${file}; done
```

4) Make table of predicted coverage for each predicted DeltaF table
```python
from make_cov_prediction_tab import *
dgFiles = [f for f in os.listdir() if 'ptDg' in f]
abundanceFile = "/hpc/hub_oudenaarden/edann/hexamers/genomes_kmers/WBcel235.kmerAbundance.csv"
tab = make_cov_prediction_tab(dgFiles, abundanceFile)
tab.to_csv('bootstrap_predictedCov.csv')
```
