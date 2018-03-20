### OBSERVED VS EXP COVERAGE
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(ggrepel)
library(pheatmap)
options(stringsAsFactors = FALSE)

pred <- read.delim('~/mnt/edann/hexamers/VAN1667prediction/predictedCoverage_avgVAN1667.txt', row.names = 1)
obs <- read.csv('~/mnt/edann/hexamers/VAN1667prediction/BS_mouse_coverage.csv', row.names = 1)

plot.obsVSexp.coverage <- function(predicted.cov, observed.cov){
  cov.df <- data.frame(template=rownames(predicted.cov), obs=observed.cov, exp=predicted.cov$exp, err=predicted.cov$err) %>% 
    mutate(lab=ifelse( exp > 20000 | obs >40000, template, '')) 
  p <- ggplot(cov.df, aes(obs, exp)) + geom_point() +  geom_text_repel(aes(label=lab),cex=4) +
    xlab('observed coverage') + ylab('predicted coverage') +
    geom_errorbar(aes(ymin=exp-err, ymax=exp+err)) +
    # ggtitle("OUD2086_distal (predicteg DeltaG on VAN1667 average)") +
    theme(axis.title = element_text(size = 20), )
    # xlim(0,2500) + ylim(0,5000)
    # geom_errorbar(aes(ymin=obs-exp-err, ymax=obs-exp+err), width=10) 
  return(p)
  }
ggsave('AvOwork/output/deltaGprediction/predictedCov_VAN1815_L2.pdf')

pred %>% mutate(lab=ifelse( exp > 20000 | obs >40000, template, '')) %>%
  ggplot(., aes(obs, exp)) + geom_point() +  # geom_text_repel(aes(label=lab),cex=2) +
  xlab('observed coverage') + ylab('obs - predicted coverage') +
  geom_errorbar(aes(ymin=exp-err, ymax=exp+err)) 
  ggtitle("OUD2086_distal (predicteg DeltaG on VAN1667 average)") +
  # xlim(0,2500) + ylim(0,5000)
  # geom_errorbar(aes(ymin=obs-exp-err, ymax=obs-exp+err), width=10) 
  ggsave('AvOwork/output/deltaGprediction/predictedCov_OUD2086_distal_residual.pdf')

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

## IS THE COVERAGE OF A TEMPLATE SEQUENCE DEPENDANT ON THE TEMPLATE ABUNDANCE?
ab <- read.csv("~/mnt/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv", 
               header=FALSE,
               col.names = c('template', 'count'))

### RSE CALCULATION
res.std.err <- function(obs,exp){ # Make sure templates are matched
  rse <- sqrt((1/length(obs))*sum((obs-exp)^2))
  return(rse)
}

pred <- pred[match(rownames(obs), rownames(pred)),]
pred.dataframe <- obs %>% mutate(exp=pred$exp, err=pred$err)
rse <- sapply(obs, function(x) res.std.err(x, pred$exp))
nreads <- colSums(obs)
df <- data.frame(sample = names(rse), n.reads=abs(nreads-mean(nreads[-1])), rse=rse)
ggplot(data = df, aes(y=rse,x=n.reads, label=sample)) + 
  theme_classic() +
  geom_point() + geom_text_repel(cex=3) +
  xlab("No. of reads - avg no. reads training set") +
  ylab("Residual Standard Error") +
  theme(axis.title = element_text(size = 18))
ggsave("~/AvOwork/output/deltaGprediction/deltaNumberOfReads.pdf")
  
## RSE for different avg
pred <- read.csv("~/mnt/edann/hexamers/zf_prediction/zf_covPrediction_avgOnNreads.csv", row.names = 1)
rse.nreads <- sapply(obs, function(sample) sapply(pred[-1], function(predDg) res.std.err(sample, predDg)))
long.rse.nreads <- melt(rse.nreads, varnames = c('training.set', 'sample'), value.name = 'RSE') %>%
  mutate(num.reads = nreads[sample], training.set=as.numeric(substr(training.set,4,100)))
ggplot(long.rse.nreads, aes(training.set, RSE, group=sample, color=num.reads)) + 
  theme_classic() +
  geom_point() + geom_line() +
  ylab('Residual Standard Error') + xlab("# of averaged samples") +
  theme(axis.title = element_text(size = 20), legend.title = element_text(size=14), legend.text = element_text(size=10))
  
ggsave("~/AvOwork/output/deltaGprediction/no.reads_in_avg_zf.pdf")
