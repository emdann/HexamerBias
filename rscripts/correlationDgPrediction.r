library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

loadPtPairs <- function(folder){
  files <- list.files(paste0(folder, "ptCounts"))
  ptPairs <- fread(paste0('mnt/edann/hexamers/rnaseq/mouse/testing/SvdB11d',no,'-MitoTrackerThird-Satellites-AdultcommonPtPairs_qual_allCells_parallel.csv'), header=TRUE, sep=',')
  colnames(ptPairs) <- c('pt', paste0('cell',colnames(ptPairs[,-1]), '.',no))
  return(ptPairs)
}
ptPairsZF <- fread(paste0('mnt/edann/hexamers/rnaseq/commonPtPairs_qual_allCells_parallel.csv'), header=TRUE, sep=',')
colnames(ptPairsZF) <- c('pt', paste0('cell',colnames(ptPairsZF[,-1]), '.zf'))

ptPairs1 <- loadPtPairs(1)
ptPairs2 <- loadPtPairs(2)
ptPairs3 <- loadPtPairs(3)
ptPairs4 <- loadPtPairs(4)
mergedPtPairs <- merge(ptPairs1, ptPairs2, by = 'pt', all=TRUE)
mergedPtPairs <- merge(mergedPtPairs, ptPairs3, by = 'pt', all=TRUE)
mergedPtPairs <- merge(mergedPtPairs, ptPairs4, by = 'pt', all=TRUE)
mergedPtPairs <- merge(mergedPtPairs, ptPairsZF, by = 'pt', all=TRUE)

# Correlation between cells
p <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) 
p1 <- p %>% mutate(GCcontPrimer=sapply(substr(ptPairs$pt, 8, 13), GCcont))
p1 %>% ggplot(.,aes(cell74,cell355, color=gccontTempl)) + geom_point() +
  scale_color_gradient2(midpoint = 0.5) +
  xlim(c(-10,-2)) + ylim(c(-10,-2)) 
  ggsave('AvOwork/output/deltaGprediction/exCor_gcTempl.pdf')

corr <- cor(mergedPtPairs[diagIx,-1],use = 'pairwise.complete.obs')
pdf('AvOwork/output/deltaGprediction/PCC_hist_gk2a_allCells.pdf')
hist(corr, breaks = 15, xlim=c(0,1), xlab='PCC', main = 'Correlation between predicted DeltaG values in single cells')
dev.off()

loadNoReads <- function(no){
  no.reads <- fread(paste0('mnt/edann/hexamers/rnaseq/mouse/SvdB11d',no,'-MitoTrackerThird-Satellites-Adult_alignedReadsPerCells.txt'), header=FALSE, sep='\t', col.names = c('no.reads','cell'))
  no.reads$cell <- paste0('cell', no.reads$cell, '.', no)
  return(no.reads)
}
no.reads.1 <- loadNoReads(1)
no.reads.2 <- loadNoReads(2)
no.reads.3 <- loadNoReads(3)
no.reads.4 <- loadNoReads(4)
no.reads.zf <- read.delim('mnt/edann/hexamers/rnaseq/gk2a-2_alignedReadsPerCells.txt', header=F, sep=' ', col.names = c('no.reads','cell'))
no.reads.zf$cell <- paste0('cell', no.reads.zf$cell, '.zf')
merged.no.reads <- rbind(no.reads.1, no.reads.2, no.reads.4, no.reads.zf)
merged.no.reads <- merged.no.reads %>% mutate(no.reads.bins = cut(merged.no.reads$no.reads, breaks = my.breaks))

my.breaks = c(seq(0,400000, by = 10000), max(merged.no.reads$no.reads))
anno.row <- data.frame(row.names =paste0('cell',no.reads.zf[match(rownames(corr), paste0('cell',no.reads$V2, '.zf')),]$V2), no.reads=no.reads[match(rownames(corr), paste0('cell',no.reads$V2),'.zf'),]$V1)

# anno = data.frame(row.names = rownames(anno.row), alignedReads=ifelse(anno.row$no.reads>40000, '> 40000', '< 40000'))

# my.breaks <- c(10000,20000,30000,max(anno.row$no.reads))
# ii <- cut(anno.row$no.reads, breaks = my.breaks, include.lowest = TRUE)
colors <- brewer.pal(3, 'RdBu')[ii]
anno <- data.frame(row.names = colnames(corr), 
                   sample = substr(colnames(corr),start=nchar(colnames(corr)), stop = nchar(colnames(corr))),
                   `n.reads(x10k)` = as.numeric(merged.no.reads$no.reads.bins[match(rownames(corr), merged.no.reads$cell)]))
p <- pheatmap(corr, show_rownames = F, show_colnames = F,
         #legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
         cellheight=2,cellwidth=2, 
         border_color = NA ,annotation_row = anno, annotation_col = anno, annotation_names_row = FALSE, annotation_names_col = FALSE, #annotation_colors = list(sample=sample,no.reads=colors),
         annotation_legend = TRUE,
         filename = 'AvOwork/output/deltaGprediction/mousezfcellsCorrDiagonal_qual_heatmap_withAnno.pdf', height = 10, width = 10)

highCov <- paste0('cell',no.reads$V2[no.reads$V1>40000])

p <- pheatmap(corr, show_rownames = F, show_colnames = F,
              legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
              cellheight=3,cellwidth=3, 
              border_color = NA,annotation_row = as.data.frame(cutree(p$tree_col, k = 5)), annotation_col = anno, annotation_names_row = FALSE, annotation_names_col = FALSE,# annotation_colors = list(no.reads=colors),
              annotation_legend = TRUE)#, fontsize = 12,

niceCl <- names(which(cutree(p$tree_col, k = 5)==1))

# Distribution of values for couple
ptOrdGcTempl <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) %>% arrange(gccontTempl)
ptDiag <- mergedPtPairs %>% filter(substr(pt,1,6)==substr(pt,8,13))
diagIx <- which(mergedPtPairs$pt %in% ptDiag$pt)

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
  y=as.numeric(ptOrdGcTempl[row,])
  par(mfrow=c(1,2), mar=c(5, 4, 4, 2) - 0.1)
  hist(y, xlab='predictedDg', breaks=15, main = ptOrdGcTempl[row,1])
  abline(v = mean(y, na.rm = T), col='red', xpd=F)
  legend('top', 'mean', lty=1, col = 'red', bty='n', xpd=T, inset = c(0,-0.05))
  qqnorm(y)
  qqline(y, xpd=F)
  test = shapiro.test(y)
  legend('topleft', paste0('SW test p-val = ', round(test$p.value, digits = 4)), cex = 0.7, bty='n')
  }

pdf('AvOwork/output/deltaGprediction/complRows_DgDistribution_qual.pdf', onefile = T, width = 10 )
sapply(complRows, rowHist, ptOrdGcTempl=ptOrdGcTempl)
dev.off()
pdf('AvOwork/output/deltaGprediction/complRows_DgDistribution_over4kreads.pdf', onefile = T, width = 10 )
sapply(complRows, rowHist, ptOrdGcTempl=ptOrdGcTempl[,c('pt',highCov, 'gccontTempl')])
dev.off()
pdf('AvOwork/output/deltaGprediction/complRows_DgDistribution_niceCl.pdf', onefile = T, width = 10 )
sapply(complRows, rowHist, ptOrdGcTempl=ptOrdGcTempl[,c('pt',niceCl, 'gccontTempl')])
dev.off()

# # Predicted deltaG
# avgDeltaG <- function(ptPairs, no.reads){
#   ptWithDg <- ptPairs %>% 
#     select(c('pt', no.reads$cell[no.reads$no.reads > 40000])) %>% 
#     mutate(predictedDg = apply(.[,-1],1, mean, na.rm=TRUE), 
#            sd = apply(.[,-1],1, sd, na.rm=TRUE))
#   predictedDg <- ptWithDg[,c('pt', 'predictedDg', 'sd')]
#   return(predictedDg)
#   }
# 
# avgDeltaG(ptPairsZF, no.reads.zf)
# write.csv(file = 'AvOwork/test_predictedDg_over40kreads.csv',x=pred)
