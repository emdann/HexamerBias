### OBSERVED VS EXP COVERAGE
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

pred <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2.CovPred.cell194.thresh.txt')
pred %>% ggplot(., aes(obs, obs - exp, label=template)) + geom_point() + geom_text(cex=2)

cells <- gsub(gsub(pattern = '.+CovPred.', replacement = '', list.files('mnt/edann/hexamers/rnaseq/predictedCov/', full.names = TRUE)), pattern = '\\..+', replacement = '')
files <- lapply(list.files('mnt/edann/hexamers/rnaseq/predictedCov/', full.names = TRUE), function(file) read.delim(file))
cors <- sapply(files, function(f) cor(f$obs, f$exp))

no.reads <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2_alignedReadsPerCells.txt', header=F, sep=' ')
no.reads$V2 <- paste0('cell',no.reads$V2)
df <- data.frame(cell=paste0('cell',cells), PCC=as.numeric(cors))
df %>% mutate(no.reads = no.reads[match(df$cell, no.reads$V2),]$V1) %>%
  ggplot(.,aes(PCC,no.reads, label=cell)) + geom_text(cex=3) +
  xlab('PCC (obs VS exp coverage)') + ylab('# aligned reads') +
  ggsave('AvOwork/output/deltaGprediction/obsexpPCCVSnoReads_noThresh.pdf')

pred <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2.CovPred.cell194.thresh.txt')
