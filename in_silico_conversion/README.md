#### IN SILICO BISULFITE CONVERSION
Scripts to compute kmer counting on in silico converted genome based on reference methylome.

## Making in silico reference genome
Script to make a reference converted genome based on deeply sequenced reference. I assume total conversion of CH sites.

1) Make binary CG reference
```
for chr in $(seq 1 19); do echo "zcat met_extraction/merged_reference_CX.cov2c.chr${chr}.gz | grep -e 'CG.' | awk '{if(\$4>\$5){met=1}else if(\$4<\$5){met=0}else{next}}{print \$1,\$2,\$3,met}' | tr ' '  '\t' > met_extraction/merged_reference_CG.met.chr${chr}"; done


```

 2)
