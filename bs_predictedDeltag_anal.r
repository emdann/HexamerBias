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

corr <- read.csv("~/mnt/edann/hexamers/zf_prediction/correlationMat_bs.csv")
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
anno <- data.frame(row.names = row.names(n.reads), num.reads=as.numeric(n.reads$num.reads))

pheatmap(corr, show_rownames = T, show_colnames = T,
         # legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
         annotation_col = anno, filename = 'AvOwork/output/deltaGprediction/BScorrelation_deltaG.pdf')

# Distribution of DeltaFs 
colnames(bs[,-1]) <- gsub(colnames(bs[,-1]), pattern = 'sorted_|_trim.', replacement = '')

long <- melt(bs, id.vars = 'TemplatePrimer_pair', variable.name = 'sample', value.name = 'predictedDg')
long %>% mutate(num.reads = as.numeric(num.reads$num.reads[num.reads$cell == as.character(long)])
ggplot(long, aes())

