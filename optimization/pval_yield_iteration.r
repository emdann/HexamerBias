suppressPackageStartupMessages(library(argparse))
library(rtracklayer)
library(purrr)
library(zoo)
library(flux)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/compare_peaks.r")

parser <- ArgumentParser()
parser$add_argument("bestEvenTrack", type="character",
                    help = "RDS file of best and even coverage predicted SCALED coverage track (the whole thing or just a sample)")
parser$add_argument("roiBed", type="character",
                    help = "BED file of regions on which to compare yields")
parser$add_argument("-t","--threads", type="integer",
                    help = "Number of cores to use")
args <- parser$parse_args()

bestEven.track.file <- args$bestEvenTrack
ROI.track.file <- args$roiBed
cores <- args$threads 

## Load track for best and even 
scaled.track <- readRDS(bestEven.track.file)

## Load track for regions of interest
roi.track <- import(ROI.track.file, format = 'BED')

## Compute p-val
score <- delta.yield.permutation(scaled.track, roi.track, threads = cores)

## Write to file
smp.name <- gsub(pattern = '.+/|\\.RDS', replacement = '', bestEven.track.file)
roi.name <- gsub(pattern = '.+/|\\.bed', replacement = '', ROI.track.file)
write(score, file = paste0('/hpc/hub_oudenaarden/edann/pred_coverage_primer_batch_D3R/evenNreads/coverage_yield/yield_pval.',smp.name,'.' ,roi.name,'.txt'), append = T)


