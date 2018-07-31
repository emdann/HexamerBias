suppressPackageStartupMessages(library(argparse))
library(rtracklayer)
library(purrr)
library(zoo)
library(flux)
library(dplyr)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/artificial_coverage/compare_peaks.r")

parser <- ArgumentParser()
parser$add_argument("tracksDir", type="character",
                    help = "path to folder containing bigWig files of tracks to compute yield")
parser$add_argument("roiBed", type="character",
                    help = "BED file of regions on which to compare yields")
# parser$add_argument("-t","--threads", type="integer",
#                     help = "Number of cores to use")
args <- parser$parse_args()

tracks.dir <- args$tracksDir
roi.track.bed <- args$roiBd

roi.track <- import(roi.track.bed, format = 'BED')

# tracks.dir <- '~/mnt/edann/pred_coverage_primer_batch_D3R/evenNreads/'

track.files <- list.files(tracks.dir, pattern = '.bw', full.names = T)

yield.df <- data.frame()
for (track.file in track.files) {
  raw.track <- import(track.file, format = 'BigWig')
  yield <- norm.scale.n.yield(raw.track, roi.track)
  yield.df <- bind_rows(yield.df, data.frame(track=gsub(track.file, pattern = '.+/', replacement = ''),
                                             yield=yield,
                                             ROI=gsub(roi.track.bed, pattern = '.+/', replacement = '')))
  }

write.csv(x = yield.df, file=paste0(tracks.dir, 'yield_', gsub(pattern = '.bed', replacement = '',as.character(yield.df$ROI[1])), '.csv'), 
          row.names = F)