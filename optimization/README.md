## Optimization of primer pool for targeted sequence enrichment

***
### Contents
* __compute_kmers_FC.py__: to compute log2(FC) between kmer composition in ROI and rest of the genome (called by wrapper ```get_kmers_ROI```)
* __compute_yield_track_dir.r__: R script to compute predicted coverage yield in regions of interest
* __FC_optimization_final.Rmd__: R notebook for E-score analysis
* __fc_optimization_functions.r__: R functions for E-score analysis (data wrangling, visualization, permutation test)
* __get_kmers_ROI.sh__: to compute for regions of interest kmer abundance, kmer abundance in random regions and log2(FC) between the two
* __pval_yield_iteration.r__: one iteration for p-val calculation of yield (called by wrapper `run_yield_pval.sh`)
* __run_combo_density.r__: scripts to compute template density for all batches with even nucleotide composition (for mouse BS-seq)
* __run_yield_pval.r__: script to compute p-value for permutation test on yield (called by wrapper `run_yield_pval.sh`)
* __run_yield_pval.sh__: wrapper to compute p-value for permutation test on yield
* __save_sample_bestVSeven_track.r__: saves sample of best VS even (random batch) regions to compute p-value on (called by wrapper `run_yield_pval.sh`)
* __which_density.Rmd__: exploratory R notebook for optimization (very initial observations)

***
### How to compute E-score
(Not tested)
1. Compute log(FC) of hexamer composition between region of interest (ROI) and whole genome **(from folder of annotation of ROI file!**)
```
cd /path/to/roi/file/
get_kmers_ROI.sh <ROI_annotation.bed> <reference_genome.fa>
```
2. Compute predicted template density for primer batch (given a set of association constants, for BS or WGS). In R:
```
source("/path/to/repo/HexamerBias/fc_optimization_functions.r")
dens.table <- density.combo(prob.vec, keqs.df, eps, nreads)
```
See R documentation for further information on input objects.

3. Compute E-score
```
roi.fc <- read.csv("/path/to/roi/file/roi.kmersFC.csv", header=F, col.names = c("template", "fc"))
fc.scores(dens.table, fc.df = roi.fc)
```

See `FC_optimization_final.Rmd` for working example.

### How to compute coverage yield over regions of interest
1. Save predicted coverage for given batch (see `FC_optimization_final.Rmd`)
2. Make predicted coverage track for batch of interest (as described in `HexamerBias/artificial_coverage/README.md`)
```
python /path/to/repo/HexamerBias/artificial_coverage/strand_specific_artificial_coverage.py <refgen.kmerAbundance.csv> <pred_coverage.csv> <random_regions.bed> /hpc/hub_oudenaarden/edann/genomes/mm10/mm10.fa --BS yes -t 15
```
3. Compute yield (giving directory containing predicted coverage track(s))
```
Rscript compute_yield_track_dir.r --threads 10 </path/to/tracks/dir/> <ROI_annotation.bed>
```

To test for significance via permutation test use `run_yield_pval.sh` (**N.B.** needs RDS file of best VS random predicted track).
