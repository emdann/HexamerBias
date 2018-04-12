#### ARTIFICIAL COVERAGE PROFILE
From the predicted coverage I want to make a igv like coverage track based on density of coverage for every hexamer (C/T).

Useful python functions are in module ```cov_from_density```.

## Comparison of artificial and real coverage
Make base resolution genome wide artificial coverage in bigWig file takes too much and comparing profiles from a
subsample would be as informative.

Take random regions of the genome (for about 1% of the mouse genome)
```
bedtools random -g /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.genome -n 500 -l 5000 -seed 42 | sort -k 1,1 -V | cut -f 1,2,3
```
Compute artificial coverage for bedfile of random regions
```
python /hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/get_artificial_cov_from_bed.py /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz predictedCoverage_avgVAN1667.txt /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/mm10.random.40.bed /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa --output bigWig
```
I can finally compared the artifical coverage track with the BigWig file of the real coverage obtained using ```bamCoverage``` from the ```deepTools``` suite.

After reading [this](https://bioconductor.org/packages/3.7/bioc/vignettes/similaRpeak/inst/doc/similaRpeak.html) I decided to use Spearman correlation as a metric for profle similarity, which assesses how well the relationship between the two profiles can be described using a monotonic function, detecting regions with possibly similar profiles even if they have different amplitudes.

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
