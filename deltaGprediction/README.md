## DNA-DNA BINDING MODEL
Scripts for implementation of model for primer binding, including
* Data preprocessing (counting genomic kmer abundance, construction of primer-template count matrix) (in bash and python)
* estimation of association constants and scaling factor epsilon (in R)
* coverage prediction (in R)

***

### Contents
* __old_modelling__: includes scripts for the first version of the model (delta G prediction)
* __Rnotebooks__: includes all notebooks for exploratory data analysis, modelling and data visualization on different datasets
* __binding_model_functions.r__: helper R functions for implementation of the model (initial data manipulation, modelling, plotting)
* __bsPrimerTemplateTab.py__: builds primer-template tab from fasta file of untrimmed reads and fasta of template sequences
* __build_pt_matrix.sh__: wrapper for construction of primer-template matrix (calls the necessary python scripts)
* __epsilon_downsampling.r__: script for epsilon estimation on downsampled bam files
* __getPrimedRegion.py__: extracts template sequences from bam files, in fasta format
* __kmersInGenome.sh__: calls kmer-counter to compute hexamer abundance in reigon of interest
* __run_cov_prediction_BS.r__: runs coverage prediction in BS-seq samples (makes K estimation on random batch and does prediction for all primer concentration)
* __run_epsilon_estimation.r__: makes chi-square estimation of epsilon from .RDS table for sample

***
