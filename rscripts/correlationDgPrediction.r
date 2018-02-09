library(dplyr)
library(data.table)
library(ggplot2)

ptPairs <- fread('mnt/edann/hexamers/rnaseq/gk2a-2.ptPairs.csv', header=TRUE, sep=',')
colnames(ptPairs) <- c('pt', paste0('cell',colnames(ptPairs[,-1])))

# Correlation between cells
p <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) 
p1 <- p %>% mutate(GCcontPrimer=sapply(substr(ptPairs$pt, 8, 13), GCcont))
p1 %>%  ggplot(.,aes(cell363,cell121, color=GCcontPrimer)) + geom_point() +
  scale_color_gradient2(midpoint = 0.5) +
  xlim(c(-10,-2)) + ylim(c(-10,-2)) +
  ggsave('AvOwork/output/deltaGprediction/exCor_gcPrimer.pdf')

corr <- cor(ptPairs[,-1],use = 'pairwise.complete.obs')
pdf('AvOwork/output/deltaGprediction/PCC_hist_gk2a.pdf')
hist(corr, breaks = 15, xlim=c(0,1), xlab='PCC', main = 'Correlation between predicted DeltaG values in single cells')
dev.off()

# Distribution of values for couple
ptOrdGcTempl <- ptPairs %>% mutate(gccontTempl=sapply(substr(ptPairs$pt, 1, 6), GCcont)) %>% arrange(gccontTempl)
ptDiag <- ptOrdGcTempl %>% filter(substr(pt,1,6)==substr(pt,8,13))
diagIx <- which(ptPairs$pt %in% ptDiag$pt)

sampleIxFromBins <- function(ptOI){
  binnedPt <- split(ptOI, ceiling(seq_along(ptOI)/100))
  smp <- sapply(binnedPt, function(bin) sample(bin, size=1))
  return(smp)
  }

boxplotPtDg <- function(ptOI, GCord=TRUE){
  if(GCord==TRUE){
    print('noice')
    ptTab <- ptOrdGcTempl
  }else{
    ptTab <- ptPairs
    print('not noice')
  }
  par(mar = c(5, 7, 4, 2))
  boxplot(t(sapply(ptTab[ptOI,-1], as.numeric)), 
          names=ptTab[ptOI]$pt,las=2, cex.names=0.4, horizontal=TRUE, cex.axis=0.7,
          xlab='predicted DeltaG', mar=c(1,2,3,10), ylim=c(-11,-2))
  } 

complRows <- which(complete.cases(ptPairs))
GCsortedRows <- complRows[order(sapply(substr(ptPairs[complRows,]$pt,1,6), GCcont))]
