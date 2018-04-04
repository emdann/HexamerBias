### AM I TAKING THE WRONG HEXAMER POSITION FOR READS MAPPING IN THE BOTTOM STRAND?
(Remember I am taking only R1 reads)
I don't think so.
- All R1 have flag 83-99 that both seem to have forward orientation in IGV. Why? Also all R2 have flag 147-163 that is said [here](https://github.com/FelixKrueger/Bismark/issues/151) to coincide to OT and OB. So it seems as if my reads are only mapping to OT and OB. Is that correct?
- The trimmed part is definetly right (checked by manual inspection of sam and fastq)

I will try and map in single-end mode all the samples in VAN1667.
```
for file in L*_R1.fastq.gz; do echo "/hpc/hub_oudenaarden/edann/bin/TrimGalore-0.4.3/trim_galore --clip_R1 9 --three_prime_clip_R1 3 --path_to_cutadapt /hpc/hub_oudenaarden/edann/venv2/bin/cutadapt -o se_mapping/ $file"; done
/hpc/hub_oudenaarden/bin/software/bismark_v0.16.3/bismark --samtools_path /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.3.1 --path_to_bowtie /hpc/hub_oudenaarden/bin/software/bowtie2-2.2.9 /hpc/hub_oudenaarden/avo/BS/mm10 $file
```
