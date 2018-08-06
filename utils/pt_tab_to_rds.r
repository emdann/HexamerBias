# sapply(c("data.table", 'reshape2', 'ggplot2', 'RColorBrewer', 'gtools', 'ggrepel'), install.packages)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Process primer-template occurrency matrix into long dataframe, saved in RDS file for quick loading in R. Includes information about kmer abundance and MFold binding energies.')
parser$add_argument("ptFile", type="character",
                    help = "File of primer-template occurrency table (.csv)")
parser$add_argument("abundanceFile", type="character",
                    help="File of genomic kmer abundance for the reference genome (.csv)")
parser$add_argument("-o", "--outputFile", type="character",
                    help="Path and name of output file")
args <- parser$parse_args()

pt.file <- args$ptFile
ab.file <- args$abundanceFile
outfile <- args$outputFile

# args <- commandArgs(trailingOnly = T)
# # Make parsing better
# pt.file <- args[1]
# ab.file <- args[2]
# outfile <- args[3]

print("Loading pt file")
pt.df <- load.pt.data(pt.file)
genome.abundance <- load.kmer.abundance(ab.file)
mfold.dg <- load.modelled.deltaG("/hpc/hub_oudenaarden/edann/hexamers/rand_hex_deltaG_ions.txt.gz")

print("Merging")
pt.all.df <- join.pt.data(pt.df$matches, pt.df$t.usage, genome.abundance, mfold.dg)

print("Saving")
saveRDS(pt.all.df %>% filter(pt>0), file = outfile)
