#!/bin/bash

smp_name=$1
ROI_bed=$2

bin_dir=/hpc/hub_oudenaarden/edann/bin/coverage_bias/optimization
R_bin=/hpc/hub_oudenaarden/edann/bin/R-3.5.1/bin

# Make sample track RDS file
${R_bin}/Rscript ${bin_dir}/save_sample_bestVSeven_track.r 500000 $smp_name

# Make a lot of iterations
for n in $(seq 1 1 1000)
  do
    echo ${R_bin}/Rscript ${bin_dir}/pval_yield_iteration.r -t 10 ${smp_name}.RDS $ROI_bed | \
      qsub -cwd -N test_iteration -pe threaded 10 -l h_rt=02:00:00
  done
