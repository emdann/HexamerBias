### Chi-square estimation of espilon for downsampling ###
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

args <- commandArgs(trailingOnly = T) 
# Make parsing better
rds.file <- args[1]

pt.df <- readRDS(rds.file)
pt.diag.df <- filter(pt.df, template==primer)
eps <- epsilon.minimize.chisq(pt.diag.df, max=200000, plot=F)

write(paste(gsub(rds.files[1], pattern='.+/|.RDS', replacement = ''), eps, sep=','), file = 'bootstrap_epsilon.txt', append = T)
