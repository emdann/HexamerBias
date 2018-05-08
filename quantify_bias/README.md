### QUANTIFICATION AND VISUALIZATION OF COVERAGE BIASES

## Bias for TSS and exons
From the coverage bigWig file obtained with ```bamCoverage``` (see artificial coverage) I make a profile of the most covered regions, using deepTools.
```
computeMatrix scale-regions \
  -R  regions_of_interest.bed \
  -S sample_coverage.bw \
  -b 3000 -a 3000 --regionBodyLength 5000 --skipZeros \
  -o sample_coverage.mat.gz
```
make the profile with deepTools:
```
plotProfile -m sample_coverage.mat.gz  \
              -out outfile.png \
              --numPlotsPerRow 1 \
              --yAxisLabel "coverage" \
              --regionsLabel samplename
```
or use functions in Rscript ```biasProfiles.r```

One liner to get average over columns (ready to plot in R)
```
zcat highcov.random.42CTCF.flank60.ratio.coverage.artCov.CTCF.mat.gz  |
  awk '{for (i=1;i<=NF;i++) if (i>=7) printf("%s ", $i); print ""}' | # Removes first 6 cols
  awk '{ for(i=1;i<=NF;i++) {total[i]+=$i ;} }
  END {
     for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;
     printf "\n" }' >  highcov.random.42CTCF.flank60.ratio.coverage.artCov.CTCF.profile.txt
```

## Comparison with non-PBAT BS-seq
The best I could find now is the data from the scWGBS protocol (Farlik et al. 2014). Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65196
