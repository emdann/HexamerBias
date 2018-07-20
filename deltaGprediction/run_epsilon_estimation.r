### Chi-square estimation of epsilon for downsampling ###
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("rdsFile", type="character",
                    help = ".RDS file of pt table")
parser$add_argument("--type", default="noBS", type="character",
                    help="Specify if sample is BS or noBS")
parser$add_argument("-o", "--outputFile", type="character",
                    help="path and name of output file")
args <- parser$parse_args()

rds.file <- args$rdsFile
type <- args$type

source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

print("Loading pt table")
pt.df <- readRDS(rds.file)
pt.diag.df <- filter(pt.df, template==primer)
if (type=="BS") {
  pt.diag.df <- pt.diag.df %>%
    filter(!grepl(template, pattern = "C"))
  }
print("Estimating espilon")
eps <- epsilon.minimize.chisq(pt.diag.df, max=200000, plot=F)

print("Writing result to file")
write(paste(gsub(rds.file, pattern='.+/|.RDS', replacement = ''), eps, sep=','), file = args$outputFile, append = T)
