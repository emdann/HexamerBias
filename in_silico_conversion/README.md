## IN SILICO BISULFITE CONVERSION
Scripts for in silico bisulfite conversion of genome based on reference methylome.

### Contents
* `convert_refgen.py`: script to build in silico converted version of reference genome based on reference methylome (see -h for details)
* `in_silico_conv_prediction.Rmd`: R notebook for prediction on coverage accounting for BS conversion
* `weighted_kmers.py`: does kmer counting in genome based on reference methylome (original script, very old, takes very long. Preferable to do the conversion with `convert_refgen.py` and then kmer counting)
