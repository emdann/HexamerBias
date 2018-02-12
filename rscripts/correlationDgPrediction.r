library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

ptPairs <- fread('mnt/edann/hexamers/rnaseq/commonPtPairs_allCells_adj.csv', header=TRUE, sep=',')
colnames(ptPairs) <- c('pt', paste0('cell',colnames(ptPairs[,-1])))

# Correlation between cells
p <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) 
p1 <- p %>% mutate(GCcontPrimer=sapply(substr(ptPairs$pt, 8, 13), GCcont))
p1 %>%  ggplot(.,aes(cell363,cell121, color=GCcontPrimer)) + geom_point() +
  scale_color_gradient2(midpoint = 0.5) +
  xlim(c(-10,-2)) + ylim(c(-10,-2)) +
  ggsave('AvOwork/output/deltaGprediction/exCor_gcPrimer.pdf')

corr <- cor(ptPairs[,-1],use = 'pairwise.complete.obs')
pdf('AvOwork/output/deltaGprediction/PCC_hist_gk2a_allCells.pdf')
hist(corr, breaks = 15, xlim=c(0,1), xlab='PCC', main = 'Correlation between predicted DeltaG values in single cells')
dev.off()

no.reads <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2_alignedReadsPerCells.txt', header=F, sep=' ')
anno.row <- data.frame(row.names =paste0('cell',no.reads[match(rownames(corr), paste0('cell',no.reads$V2)),]$V2), no.reads=no.reads[match(rownames(corr), paste0('cell',no.reads$V2)),]$V1)
anno = data.frame(row.names = rownames(anno.row), alignedReads=ifelse(anno.row$no.reads<30000, '< 30000', '> 30000'))
my.breaks <- c(10000,20000,30000,max(anno.row$no.reads))
ii <- cut(anno.row$no.reads, breaks = my.breaks, include.lowest = TRUE)
colors <- brewer.pal(3, 'RdBu')[ii]
p <- pheatmap(corr, show_rownames = F, show_colnames = F,
         legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
         cellheight=2.5,cellwidth=2.5, 
         border_color = NA,annotation_row = anno, annotation_col = anno, annotation_names_row = FALSE, annotation_names_col = FALSE,# annotation_colors = list(no.reads=colors),
         annotation_legend = TRUE) #, fontsize = 12,
         filename = 'AvOwork/output/deltaGprediction/cellsCorr_heatmap_withAnno.pdf', height = 10, width = 10)

corClust <- names(which(cutree(p$tree_col, k=2)==1))

# Distribution of values for couple
ptOrdGcTempl <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) %>% arrange(gccontTempl)
ptDiag <- ptOrdGcTempl %>% filter(substr(pt,1,6)==substr(pt,8,13))
diagIx <- which(ptOrdGcTempl$pt %in% ptDiag$pt)

complRows <- which(complete.cases(ptOrdGcTempl))
GCsortedRows <- complRows[order(sapply(substr(ptPairs[complRows,]$pt,1,6), GCcont))]

sampleIxFromBins <- function(ptOI){
  binnedPt <- split(ptOI, ceiling(seq_along(ptOI)/100))
  smp <- sapply(binnedPt, function(bin) sample(bin, size=1))
  return(smp)
  }

boxplotPtDg <- function(ptOI, ptTab=ptOrdGcTempl, GCord=TRUE){
  if(GCord==TRUE){
    # print('noice')
    ptTab <- ptTab
  }else{
    ptTab <- ptPairs
    # print('not noice')
  }
  ## make color scale
  values=ptTab[ptOI,]$gccontTempl
  ii <- cut(values, breaks = seq(0, 1, len = 100), include.lowest = TRUE)
  colors <- colorRampPalette(c('lightblue','white','firebrick'))(99)[ii]
  par(mar = c(5, 7, 4, 2), xpd=TRUE)
  boxplot(t(sapply(ptTab[ptOI,-1], as.numeric)), 
          names=ptTab[ptOI,]$pt,las=2, cex.names=0.4, horizontal=TRUE, cex.axis=0.7, col=colors, 
          xlab='predicted DeltaG', mar=c(1,2,3,10), ylim=c(-11,-1))
  legend("topleft", inset=c(0.4,-0.05), legend=c("> AT","> GC"), fill = c('lightblue', 'firebrick'), horiz = TRUE, bty = 'n')
} 

rowHist <- function(row, ptOrdGcTempl){
  y=as.numeric(ptOrdGcTempl[row,-ncol(ptOrdGcTempl)])
  par(mfrow=c(1,2), mar=c(5, 4, 4, 2) - 0.1)
  hist(y, xlab='predictedDg', breaks=15, main = ptOrdGcTempl[row,1], xlim=c(-1,-12))
  abline(v = mean(y, na.rm = T), col='red', xpd=F)
  legend('top', 'mean', lty=1, col = 'red', bty='n', xpd=T, inset = c(0,-0.05))
  qqnorm(y)
  qqline(y, xpd=F)
  test = shapiro.test(y)
  legend('topleft', paste0('SW test p-val = ', round(test$p.value, digits = 4)), cex = 0.7, bty='n')
  }

pdf('AvOwork/output/deltaGprediction/complRows_DgDistribution.pdf', onefile = T, width = 10 )
sapply(complRows, rowHist, ptOrdGcTempl=ptOrdGcTempl)
dev.off()
pdf('AvOwork/output/deltaGprediction/complRows_DgDistribution_goodCluster.pdf', onefile = T, width = 10 )
sapply(complRows, rowHist, ptOrdGcTempl=ptOrdGcTempl[,c('pt',corClust)])
dev.off()

# Predicted deltaG
ptWithDg <- ptOrdGcTempl %>% mutate(predictedDg = apply(ptOrdGcTempl[,rownames(anno.row>30000)],1, median, na.rm=TRUE)) 
write.csv(file = 'AvOwork/predictedDg_over30kreads.csv',x=ptWithDg[c('pt', 'predictedDg')])
ptWsd <- ptWithDg %>% mutate(sd=apply(ptOrdGcTempl[,rownames(anno.row>30000)],1, sd, na.rm=TRUE)) %>% 
  ggplot(., aes(sd)) + geom_histogram()
