## PREDICTED COVERAGE PROFILES FROM PRIMER BINDING
From the template usage I want to make a IGV like coverage track based on density of coverage for every hexamer (where density is the ratio between the template coverage and the template genomic abundance).

### Contents
* `artCov_primer_variation.py`: makes predicted coverage tracks for progressive changes of random primer compositions
* `compare_peaks.r`: helper R functions to load, plot and statistical analysis of predicted VS experimental profiles
* `cov_from_density.py`: helper python functions for computation of predicted profiles
* `density_EDA.Rmd`: notebook of exploratory data analysis
* `genome_wide_artificialcov.py`: computes genome wide predicted coverage (Never run until the end, probably takes years)
* `predicted_tracks_final.Rmd`: R notebook with comparison between predicted and experimental tracks in many samples.
* `strand_specific_artificial_coverage.py`: main function to compute predicted profiles in regions of interest defined in a bed file

***

### How to make a predicted coverage profile
1. Take random regions of the genome of choice (doing whole genome will take very very long):
```
bedtools random -g /path/to/genome/file/${mygenome}.genome -n $howmanyregions -l $lengthregions -seed 42 | sort -k 1,1 -V | cut -f 1,2,3
```
The script can compute the predicted track for multiple regions in parallel, so better to have more regions (higher n) of smaller length (small l).
__CHECK FOR AND REMOVE OVERLAPS BETWEEN REGIONS!__ Otherwise the script will crash at the very end because the bigWig file has to be ordered. For example, if l=3000:
```
cat $bedfile | bedtools merge | awk '$3-$2==3000' > bedfile.noOvs.bed
```
2. Build predicted coverage track for regions:
```
python strand_specific_artificial_coverage.py ${genome_kmer_abundance}.csv ${kmer_coverage}.csv bedfile.noOvs.bed ${reference_genome}.fasta   /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa --output bigWig
```
 Give it a lot of threads to use! Run `python strand_specific_artificial_coverage.py --help` to check further options and file formats.

The output predicted track in bigWig format you can upload in R for visualization and comparison with other (experimental or predicted) coverage tracks.

<!-- #### Computing average coverage over coding regions
Using the deepTools suite:
```
computeMatrix scale-regions \
  -R genes.bed \
  -S ${pred_cov_track}.bw \
  -b 3000 -a 3000 --regionBodyLength 5000 --skipZeros \
  -o ${pred_cov_track}_genes.mat.gz
``` -->
