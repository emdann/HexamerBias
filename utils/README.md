## USEFUL SCRIPTS

### Contents
* `array_job_wrapper.sh`: wrapper script to send an array of jobs to the cluster
* `count_pt_pairs.py`: script to count the total number of primer-template pairs in the co-occurrency matrix (`.ptCounts.qualFilt.csv`). Useful for example to check how many reads passed the quality control filtering.
* `deduplicateCELseqSamfile.py`: removes PCR duplicates in CELseq2 libraries (filters out reads with the same position, umi, cell, primer sequence)
* `find_genomecov.sh` **(!!!)**: computes percentage of genome with coverage=0 calling bedtools (for extrapolation of epsilon)
* `freeEnergy.py` **(!!!)**: functions to compute free energies for binding of complementary sequences using NN model. Includes specification for ionic conditions.
* `index_bam.sh`: wrapper to sort and index bam files with Samtools
* `kmerCounter.py` **(!!!)**: script to count kmer abundance in fasta file (default kmer length: 6 bps). Runs on multiple cores (for example if you are counting in a whole genome)
* `mergeCov2c.py`: merges cytosine reports from different chromosomes in one file (output of `bismark_methylation_extractor --cytosine_report --split_by_chromosome`)
* `nuc_tss.py`: computes avg base composition for n sequences of the same length from fasta file
* `pt_tab_to_rds.r` **(!!!)**: converts primer-template matrix from `.ptCounts.qualFilt.csv` to long dataframe R object file `.RDS`, for quick data import
* `send_array_job,sh`: called by wrapper `array_job_wrapper.sh`
* `splitByChr.py`: splits cytosine reports (cov2c) by chromosomes
* `subsmpBam.py`: makes random subsample of bam file to a certain fraction of reads
* `subsmpCsv.py`: makes random subsample of rows froma a big csv file
