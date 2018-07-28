suppressPackageStartupMessages(library(argparse))
library(rtracklayer)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/compare_peaks.r")

parser <- ArgumentParser()
parser$add_argument("sampleSize", type="integer",
                    help = "Size of sample")
parser$add_argument("output", type="character",
                    help = "output file name")
args <- parser$parse_args()

sample.size <- args$sampleSize
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

smp <- sample(scaled.track, sample.size)
saveRDS(smp, file = paste0('/hpc/hub_oudenaarden/edann/pred_coverage_primer_batch_D3R/evenNreads/coverage_yield/', out.name, '.RDS'))