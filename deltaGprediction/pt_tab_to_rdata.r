install.packages("data.table", 'reshape2', 'ggplot2', 'RColorBrewer', 'gtools', 'ggrepel')
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

args <- commandArgs(trailingOnly = T) 
# Make parsing better
pt.file <- args[1]
ab.file <- args[2]
outfile <- args[3]

print("Loading pt file")
pt.df <- load.pt.data(pt.file)
genome.abundance <- load.kmer.abundance(ab.file)
mfold.dg <- load.modelled.deltaG("/hpc/hub_oudenaarden/edann/hexamers/rand_hex_deltaG_ions.txt.gz")

print("Merging")
pt.all.df <- join.pt.data(pt.df$matches, pt.df$t.usage, genome.abundance, mfold.dg)

print("Saving")
save(pt.all.df, file = outfile)