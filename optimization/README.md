## Optimization of primer pool to direct coverage

***
## Contents
* __DEimplementation.py__: to run differential evolution optimization
* __get_kmers_ROI.sh__: to compute for regions of interest kmer abundance, kmer abundance in random regions and log2(FC) between the two
* __compute_kmers_FC.py__: to compute log2(FC) between ROI and random kmers (called by ```get_kmers_ROI```)
* __primerProbability.py__: functions related to optimisation and computing probabilities and from matrix of nucleotides
* __tests_for_DE.py__: functions to test performance of optimisation algorithm
* __DEperformance.Rmd__: R notebook for analysis of DE optimisation output
* __shiny_DE__: shiny app for visualisation of DE optimisation output
