suppressPackageStartupMessages(library(argparse))
library(rtracklayer)
library(purrr)
library(zoo)
library(flux)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/compare_peaks.r")

parser <- ArgumentParser()
parser$add_argument("output", type="character",
                    help = "output file name")
args <- parser$parse_args()


out.name <- args$output

## Load track for best and even 
load('/hpc/hub_oudenaarden/edann/bestVSeven_track.RData', verbose=T)

## Scale to remove values under zero
scaled.track <- best.even.track.all
score.cols <- colnames(values(scaled.track)[sapply(values(scaled.track), is.numeric)])
min.score <- min(as.matrix(scaled.track@elementMetadata[score.cols]))
for (col in score.cols) {
  scaled.track@elementMetadata[col][[1]] <- scaled.track@elementMetadata[col][[1]] + abs(min.score)
}

smp <- sample(scaled.track, 1000000)
saveRDS(smp, file = paste0('/hpc/hub_oudenaarden/edann/pred_coverage_primer_batch_D3R/evenNreads/coverage_yield/', out.name, '.RDS'))