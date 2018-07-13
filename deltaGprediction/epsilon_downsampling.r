### Chi-square estimation of espilon for downsampling ###
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

args <- commandArgs(trailingOnly = T) 
# Make parsing better
rds.file <- args[1]

print("Loading pt table")
pt.df <- readRDS(rds.file)
pt.diag.df <- filter(pt.df, template==primer)

print("Estimating espilon")
eps <- epsilon.minimize.chisq(pt.diag.df, max=200000, plot=F)

print("Writing result to file")
write(paste(gsub(rds.file, pattern='.+/|.RDS', replacement = ''), eps, sep=','), file = 'bootstrap_epsilon.txt', append = T)
