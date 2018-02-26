### OBSERVED VS EXP COVERAGE
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

pred <- read.delim('mnt/edann/hexamers/rnaseq/mouse/predictedCov/SvdB11d1-MitoTrackerThird-Satellites-Adult.CovPred.10.thresh1.qual.txt')
pred %>% ggplot(., aes(obs, obs-exp, label=template)) + geom_point() + geom_text(cex=2) +
  xlab('observed coverage') + ylab('predicted coverage') 
  ggsave('AvOwork/output/deltaGprediction/predictedCov_avg_cell30_notbadcor.pdf')

cell2cell.files <- list.files('mnt/edann/hexamers/rnaseq/predictedCov/',pattern =  'gk2a-2.CovPred.+qual', full.names = TRUE)
cell2avg.files <- list.files('mnt/edann/hexamers/rnaseq/zebrafish/predictedCov/',pattern =  'gk2a-2.40k.', full.names = TRUE)

makeObsExpCorrelation <- function(files){
  cells <- gsub(gsub(pattern = '.+CovPred.', replacement = '', files), pattern = '\\..+', replacement = '')
  read.files <- lapply(files, function(file) read.delim(file))
  cors <- sapply(read.files, function(f) cor(f$obs, f$exp))
  return(data.frame(cell=paste0('cell',cells), PCC=as.numeric(cors)))
  }

no.reads <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2_alignedReadsPerCells.txt', header=F, sep=' ')
no.reads$V2 <- paste0('cell',no.reads$V2)
df <- data.frame(cell=paste0('cell',cells), PCC=as.numeric(cors))

plotPcc <- function(df){
  df %>% mutate(no.reads = no.reads[match(df$cell, no.reads$V2),]$V1) %>%
    ggplot(.,aes(PCC,no.reads, label=cell)) + geom_text(cex=3) +
    xlab('PCC (obs VS exp coverage)') + ylab('# aligned reads') 
#  ggsave('AvOwork/output/deltaGprediction/obsexpPCCVSnoReads_40kavg_qual.pdf')
}

resMatrix <- function(files){
  cells <- gsub(gsub(pattern = '.+CovPred.', replacement = '', files), pattern = '\\..+', replacement = '')
  read.files <- lapply(files, function(file) read.delim(file))
  residual.mat = do.call(cbind, lapply(read.files, function(f) mutate(f, res=obs-exp)$res))
  colnames(residual.mat) <- paste0('cell',cells)
  rownames(residual.mat) <- read.files[[1]]$template
  return(residual.mat)
  }
  
plotResMat <- function(resmat){
p <- pheatmap(resmat[,no.reads$V2[no.reads$V1>40000]], breaks = seq(-100,100, length.out = 100), 
                annotation_col = data.frame(row.names = no.reads$V2, no.reads=no.reads$V1),
                show_rownames = F, show_colnames = T, fontsize_col = 6)
  return(p)
  }

clust <- cutree(celltree$tree_row, k = 6)

pdf('AvOwork/output/deltaGprediction/cellHex_residuals_avg_over40kcells_annoReads.pdf', width = 10, height = 10)
plotResMat(avg.resmat)
dev.off()

pdf('AvOwork/output/deltaGprediction/cellHex_residuals_cell_over40kcells_annoReads.pdf', width = 10, height = 10)
plotResMat(cell.resmat)
dev.off()

### AVG RESIDUAL AND AVG OBS
read.files <- lapply(cell2avg.files, function(file) read.delim(file))

avgDf <- data.frame(avgObs = apply(sapply(read.files, function(x) x$exp), 1, mean), avgRes=apply(avg.resmat,1,mean))
cellDf <- data.frame(avgObs = apply(sapply(read.files, function(x) x$exp), 1, mean), avgRes=apply(cell.resmat,1,mean))

avgDf %>% mutate(template =rownames(avgDf)) %>% mutate(lab=ifelse(abs(avgRes) > 20, template, '')) %>%
  ggplot(., aes(avgObs, avgRes, label=lab)) + geom_point() + geom_text(cex=2, nudge_y = 0.5) +
  geom_hline(yintercept = 0, color='red') +
  ggtitle('Avg DeltaG') + xlab('mean(observed coverage)') + ylab('mean(observed - expected coverage)') 
  ggsave('AvOwork/output/deltaGprediction/avgResidualplot_avg_zf.pdf')

cellDf %>% mutate(template =rownames(cellDf)) %>% mutate(lab=ifelse(abs(avgRes) > 20, template, '')) %>%
ggplot(., aes(avgObs, avgRes, label=lab)) + geom_point() + geom_text(cex=2, nudge_y = 0.5) +
  geom_hline(yintercept = 0, color='red') +
  ggtitle('cell DeltaG')+ xlab('mean(observed coverage)') + ylab('mean(observed - expected coverage)') 
  ggsave('AvOwork/output/deltaGprediction/avgResidualplot_cell_zf.pdf')

