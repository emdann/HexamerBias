### QUANTIFICATION AND VISUALIZATION OF COVERAGE BIASES

## Bias for TSS and exons
From the coverage bigWig file obtained with ```bamCoverage``` (see artificial coverage) I make a profile of the most covered regions, using deepTools.
```
computeMatrix scale-regions \
  -R  /hpc/hub_oudenaarden/edann/hexamers/annotations_bed/RefSeq_genes_mm10-2.bed \ # separate multiple files with spaces
  -S /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/artificial_coverage/VAN1667.cov.bw \
  -b 3000 -a 3000 \
  --regionBodyLength 5000 \
  --skipZeros -o /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/artificial_coverage/matrixVAN1667.gz \
  --outFileNameMatrix /hpc/hub_oudenaarden/edann/hexamers/VAN1667prediction/artificial_coverage/matrixVAN1667.tab \
```
