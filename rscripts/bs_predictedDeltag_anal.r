### Predicted Delta F on bisulfite seq data
library(dplyr)
library(ggplot2)
library(data.table)

# bs <- fread("~/mnt/edann/hexamers/zf_prediction/BS.commonPtPairs.csv", sep=',', header=TRUE)
load('~/BS.commonPtPairs.RData')

## Correlation heatmaps 
numreads <- read.delim("~/mnt/edann/hexamers/VAN1667prediction/VAN1667.numReads.txt", header=FALSE, col.names=c('numreads', 'cell'), nrows = 9)
corr <- read.csv('~/mnt/edann/hexamers/VAN1667prediction/corrMatrixDg.csv', row.names = 1)
pheatmap(corr, breaks=(seq(0,1,0.01)), annotation_row = numreads[-2], filename = 'AvOwork/output/deltaGprediction/VAN1667_correlation_dG.pdf')

corr <- read.csv("~/mnt/edann/hexamers/zf_prediction/correlationMat_bs.csv", row.names=1)
colnames(corr) <- gsub(gsub(colnames(corr), pattern = 'sorted_', replacement = ''), pattern='_trim.', replacement = '')
rownames(corr) <- gsub(gsub(rownames(corr), pattern = 'sorted_', replacement = ''), pattern='_trim.', replacement = '')


mm2 <- read.delim("~/mnt/edann/hexamers/VAN1667prediction/VAN1667.numReads.txt", header=FALSE, col.names = c('num.reads', 'cell'))
zf <- read.delim("~/mnt/edann/hexamers/zf_prediction/zebrafishBs.numReads.txt", header=FALSE, col.names = c('num.reads', 'cell'))
num.reads <- rbind(zf,mm2)[-26,]
rownames(num.reads)<- num.reads$cell
rownames(num.reads) <- gsub(rownames(num.reads), pattern='_trim.|_bismark.+', replacement = '')
rownames(num.reads)[grep(rownames(num.reads), pattern = 'R.', invert = TRUE)] <- paste0(rownames(num.reads)[grep(rownames(num.reads), pattern = 'R.', invert = TRUE)], '_R1')

rownames(num.reads) <- gsub(rownames(num.reads), pattern='-', replacement = '.')
n.reads <- num.reads[match( colnames(corr), rownames(num.reads)),] 
anno <- data.frame(row.names = row.names(n.reads), 
                   num.reads=as.numeric(n.reads$num.reads), 
                    organism=ifelse(sapply(as.character(rownames(corr)), function(x) grepl(x, pattern="^L.+")), 'mouse', 'zebrafish'))

pheatmap(corr, show_rownames = T, show_colnames = T,
         # legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
         annotation_col = anno, 
         filename = 'AvOwork/output/deltaGprediction/BScorrelation_deltaG.pdf')

# Distribution of DeltaFs 
colnames(bs) <- c("TemplatePrimer_pair", gsub(colnames(bs[, -"TemplatePrimer_pair"]), pattern = 'sorted_|_trim.', replacement = ''))

num.reads$cell <- gsub(num.reads$cell, pattern = 'trim._|_bismar.+', replacement = '')
num.reads$cell[grep(num.reads$cell, pattern = 'R.', invert = TRUE)] <- paste0(num.reads$cell[grep(num.reads$cell, pattern = 'R.', invert = TRUE)], '_R1')


long <- melt(bs, id.vars = 'TemplatePrimer_pair', variable.name = 'sample', value.name = 'predictedDg')
long.bs <- long %>% mutate(num.reads=num.reads$num.reads[sample],
                           organism=ifelse(grepl(sample, pattern = '^L'), 'mouse', 'zf'))
p <- ggplot(long.bs, aes(y=predictedDg, x=num.reads)) + 
  geom_boxplot(aes(group=cut_width(num.reads, 10000)), outlier.alpha = 0.1)

### Distribution of deltaG across samples

# Read number and deltag distribution
ggplot(long.bs, aes(as.numeric(num.reads), predictedDg)) +
  geom_boxplot(aes(group=sample, color=organism), alpha=0.2, width=100000, outlier.alpha = 0.2) +
  xlab('No. of reads') + ylab("Predicted DeltaG") 
  ggsave("~/AvOwork/output/deltaGprediction/deltaGdistributionVSnumreads_bs.pdf")


smp <- sample_n(bs, size = 2000) # Random sampling of pt pairs

bs.diag <- take.diagonal(bs) # Take diagonal pt

# Make boxplots for pt pairs sorted by GC content
allsmp <- smp %>% is.na %>% rowSums()<8 # Take pt pairs that have values for both organisms
compl.smp <- smp[allsmp,]

pdf("~/AvOwork/output/deltaGprediction/ptPairsDg_organism_bs.pdf", onefile=TRUE, width = 10, height = 10)
first=1
while (first < nrow(compl.smp)) {
  long.smp <- melt(compl.smp[first:(first+35),], id.vars = 'TemplatePrimer_pair', variable.name = 'sample', value.name = 'predictedDg') %>% 
    mutate(num.reads = num.reads$num.reads[sample], 
           organism = ifelse(grepl(sample, pattern = '^L'), 'mouse', 'zf'),
           GCcontent = sapply(substr(TemplatePrimer_pair, 1,6), GCcont)) %>%
    arrange(desc(GCcontent))
  
  p <- ggplot(long.smp, aes(organism, predictedDg)) +
    facet_wrap( ~ TemplatePrimer_pair, nrow=6, ncol=6) +
    theme(strip.text = element_text(size=8)) +
    geom_boxplot(aes(fill=organism), alpha=0.2,outlier.alpha = 0.2, varwidth = TRUE) +
    ylab("Predicted DeltaG")
  print(p)
  first = first + 36
  }
dev.off()

ggsave("~/AvOwork/output/deltaGprediction/deltaGdistributionVSnumreads_bs.pdf")
