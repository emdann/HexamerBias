### AM I TAKING THE WRONG HEXAMER POSITION FOR READS MAPPING IN THE BOTTOM STRAND?
(Remember I am taking only R1 reads)
I don't think so.
- All R1 have flag 83-99 that both seem to have forward orientation in IGV. Why? Also all R2 have flag 147-163 that is said [here](https://github.com/FelixKrueger/Bismark/issues/151) to coincide to OT and OB. So it seems as if my reads are only mapping to OT and OB. Is that correct?
- The trimmed part is definetly right (checked by manual inspection of sam and fastq)


### STRAND SPECIFC PREDICTED COVERAGE
I will try and map in single-end mode all the samples in VAN1667.
```
for file in L*_R1.fastq.gz; do echo "/hpc/hub_oudenaarden/edann/bin/TrimGalore-0.4.3/trim_galore --clip_R1 9 --three_prime_clip_R1 3 --path_to_cutadapt /hpc/hub_oudenaarden/edann/venv2/bin/cutadapt -o se_mapping/ $file"; done
echo "/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3/bismark --samtools_path /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 --path_to_bowtie /hpc/hub_oudenaarden/bin/software/bowtie2-2.2.9 /hpc/hub_oudenaarden/avo/BS/mm10 $file" | qsub -cwd -N map_L1_R2 -pe threaded 10 -l h_rt=24:00:00 -l h_vmem=50G -l h_cpu=1:00:00
echo "/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3/deduplicate_bismark --samtools_path /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 -s --bam ${file}" | qsubl -N dedup

/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3/bismark_methylation_extractor --samtools_path /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 -s --gzip --report --multicore 8 --cytosine_report  --genome_folder ${refgen_dir} $file
```
Merge bam files
```
/hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools cat -o VAN1667_se.bam L*_R1_trimmed_bismark_bt2.deduplicated.bam
```

I do all the Predicted Dg computation

### Strand specific artificial coverage
To do the read extension I can't parallelize over chunks of sequences, so I better parallelize on the bed entries.

1) Make bed file of many (short) random regions of the genome

2) Run file
```
echo "source /hpc/hub_oudenaarden/edann/venv2/bin/activate; python /hpc/hub_oudenaarden/edann/bin/coverage_bias/strand_specific_artificial_coverage.py /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz /hpc/hub_oudenaarden/edann/hexamers/strand_specific/VAN1667_se_predictedcov.csv artificial_coverage/mm10.random.42.bed /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa -t 10" | qsubl -N artificial_cov -pe threaded 10
```
