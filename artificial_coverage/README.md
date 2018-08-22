## PREDICTED COVERAGE PROFILES FROM PRIMER BINDING
From the template usage I want to make a IGV like coverage track based on density of coverage for every hexamer $\frac{C}{T})$.

## Contents
* __artCov_primer_variation.py__: makes predicted coverage tracks for progressive changes of random primer compositions
* __compare_peaks.r__: helper functions to load, plot and statistical analysis of predicted VS experimental profiles in R
* __cov_from_density.py__: helper functions for profile computation
* __density_EDA.Rmd__: notebook of exploratory data analysis
* __genome_wide_artificialcov.py__: computes genome wide predicted coverage (probably takes years)
* __strand_specific_artificial_coverage.py__: main function to compute predicted profiles in regions of interest defined in a bed file

## Comparison of artificial and real coverage
Make base resolution genome wide artificial coverage in bigWig file takes too much and comparing profiles from a
subsample would be as informative.

Take random regions of the genome (for about 1% of the mouse genome)
```
bedtools random -g /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.genome -n 500 -l 5000 -seed 42 | sort -k 1,1 -V | cut -f 1,2,3
```
CHECK FOR OVERLAPS!!! (Otherwise the order in saving the wig is messed up)
```
cat $bedfile | bedtools merge | awk '$3-$2==3000' > bedfile.noOvs.bed
```

Compute artificial coverage for bedfile of random regions
```
python /hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/get_artificial_cov_from_bed.py /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz predictedCoverage_avgVAN1667.txt /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.random.40.bed /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa --output bigWig
```
I can finally compared the artifical coverage track with the BigWig file of the real coverage obtained using ```bamCoverage``` from the ```deepTools``` suite.

After reading [this](https://bioconductor.org/packages/3.7/bioc/vignettes/similaRpeak/inst/doc/similaRpeak.html) I decided to use Spearman correlation as a metric for profle similarity, which assesses how well the relationship between the two profiles can be described using a monotonic function, detecting regions with possibly similar profiles even if they have different amplitudes.

### Sampling regions with high coverage
```
source /hpc/hub_oudenaarden/edann/venv2/bin/activate; bamCoverage -b /hpc/hub_oudenaarden/edann/crypts_bs/VAN1667/se_mapping/VAN1667_se.bam -o VAN1667_se.bw -p 6 -bs 5000
bigWigToWig VAN1667.cov.bw VAN1667.cov.wig
cat VAN1667_se.wig | cut -f 4 | sort -n | uniq -c > VAN1667_se.cov.hist.txt

```

Find experimental coverage in regions of interest
```
cat myregions.bed |
  bedtools makewindows -b stdin -w 50 > myregions.50w.bed
  samtools bedcov myregions.50w.bed /hpc/hub_oudenaarden/edann/crypts_bs/VAN1667/se_mapping/VAN1667_se.srt.bam
```

### More random things

Make 100 bps bins in the bedfile used for artificial coverage
```
bedtools makewindows -b mm10.random.40.bed -w 100 > mm10.random.40.bins.bed
```

Put both (or many more) profiles on the same file
```
multiBigwigSummary BED-file --bwfiles sorted_L1_trim1_R1_bismark_bt2_pe.cov.bw mm10.random.40.artCov.bw -o test_multiBigWig.npz --BED mm10.random.40.bins.bed --outRawCounts test_multiBigWig.tab
```

## Does primer probability proportional to template abundance give even coverage?
I compute predicted coverage with primer probability proportional to template abundance with the ```get_proportional_coverage``` function in the ```primerProbability``` script.

## Find the baseline
1) Make bedgraph from bigWig files of experimental and artificial coverage
```
bigWigToWig -chrom=chr1 artificial_coverage/mm10.random.42.3000smp.artCov.bw artificial_coverage/mm10.random.42.3000smp.artCov.wig
cat artificial_coverage/mm10.random.42.1000smp.wig | tail -n+2 | awk '{print "chr1\t"$1"\t"$1+1"\t"$2}' > artificial_coverage/mm10.random.42.1000smp.bedGraph
```
2) Intersect regions covered or not covered regions in the experimental data
```
cat VAN1667_se.random42.srt.chr1.wig | awk '$4==0' | bedtools intersect -a stdin -b artificial_coverage/mm10.random.42.1000smp.bedGraph -wb | cut -f 8
cat VAN1667_se.random42.srt.chr1.wig | awk '$4>5' | bedtools intersect -a stdin -b artificial_coverage/mm10.random.42.1000smp.bedGraph -wb | cut -f 8
```
