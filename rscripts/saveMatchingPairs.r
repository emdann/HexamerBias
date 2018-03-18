library(data.table)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(hexbin)


load.common.pt.file <- function(commonPtFile, save = FALSE){
  sampleName <- gsub(gsub(commonPtFile, pattern = '.commonPtPairs.csv', replacement = ""), pattern = ".+/", replacement = '')
  ptPairs <- fread(commonPtFile, header=TRUE, sep=',')
  colnames(ptPairs) <- c('TemplatePrimer_pair', paste0(sampleName,'.cell', colnames(ptPairs)[-1]))
  if (save) {save(ptPairs, file = paste0(sampleName, '.commonPtPairs.RData'))}
  return(ptPairs)
}
  
load.numReads.file <- function(numReadsFile, save = FALSE){
  sampleName <- gsub(gsub(numReadsFile, pattern = '.numReads.txt', replacement = ""), pattern = ".+/", replacement = '')
  num.reads <- read.delim(numReadsFile, header=TRUE, sep='\t', col.names = c( 'cell', 'num.reads')) %>%
    mutate(cell=paste0(sampleName,'.cell',cell)) %>% mutate(cell=gsub(cell, pattern='-', replacement = '.'))
  if (save) {save(num.reads, file = paste0(sampleName, '.numReads.RData'))}
  return(num.reads)
} 

make.long.commonPt <- function(ptPairs, num.reads, file=NULL) {
  ## Possibly subset the pairs you need
  data <- ptPairs %>% 
    melt(id.vars = 'TemplatePrimer_pair', variable.name = 'cell', value.name = 'deltaG') %>%
    mutate(gc_content = sapply(substr(TemplatePrimer_pair,1,6), GCcont)) 
  data <- data %>% mutate(num_reads = num.reads$num.reads[match(data$cell, num.reads$cell)])
  if (!is.null(file)) {save(data, file = file)}
  return(data)
}

take.diagonal <- function(ptPairs){
  ptPairs <- ptPairs %>% filter(substr(TemplatePrimer_pair,1,6)==substr(TemplatePrimer_pair,8,13))
  return(ptPairs)
  }

logPtPairs <- data.frame(TemplatePrimer_pair=ptPairs$TemplatePrimer_pair, log(ptPairs[,-1]))

make.cell.corr <- function(ptPairs){ 
  corr <- cor(ptPairs[,-1],use = 'pairwise.complete.obs') 
  return(corr)
  }

colors <- brewer.pal(3, 'RdBu')[ii]

### Work on RNA seq data
mm1 <- load.common.pt.file("~/mnt/edann/hexamers/rnaseq/mouse/testing/SvdB11d1-MitoTrackerThird-Satellites-Adult.commonPtPairs.csv")
mm2 <- load.common.pt.file("~/mnt/edann/hexamers/rnaseq/mouse/SvB11d2/SvdB11d2-MitoTrackerThird-Satellites-Adult.commonPtPairs.csv")
zf <- load.common.pt.file("~/mnt/edann/hexamers/rnaseq/zebrafish/gk2a-2.commonPtPairs.csv")

mm1 <- load.numReads.file("~/mnt/edann/hexamers/rnaseq/mouse/testing/SvdB11d1-MitoTrackerThird-Satellites-Adult.commonPtPairs.csv")

mm1 <- load.numReads.file("~/mnt/edann/hexamers/rnaseq/mouse/testing/SvdB11d1-MitoTrackerThird-Satellites-Adult.numReads.txt")
mm2 <- read.delim("~/mnt/edann/hexamers/VAN1667prediction/VAN1667.numReads.txt", header=FALSE, col.names = c('num.reads', 'cell'))
zf <- read.delim("~/mnt/edann/hexamers/zf_prediction/zebrafishBs.numReads.txt", header=FALSE, col.names = c('num.reads', 'cell'))

n.reads <- num.reads$num.reads[match( colnames(corr), num.reads$cell)]
anno <- data.frame( `num.reads(x10k)` = as.numeric(cut((n.reads), breaks = c(seq(10000,100000, by=10000), max(n.reads))))
                  sample = gsub(colnames(corr), pattern = '.cell.+', replacement = ''),
                  row.names = colnames(corr)) 
# anno <- anno %>% mutate(organism = as.factor(ifelse(anno$sample=='gk2a.2', 'zf', 'mouse')))
pheatmap(corr, show_rownames = F, show_colnames = F,
              # legend_breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, max(corr)), legend_labels=c(0.5, 0.6, 0.7, 0.8, 0.9, 'PCC'),
              cellheight=0.6, cellwidth=0.6) ,
              border_color = NA ,annotation_row = anno,  annotation_col = anno,
              annotation_names_row = FALSE, annotation_names_col = FALSE, #annotation_colors = list(sample=sample,no.reads=colors),
              annotation_legend = TRUE,
              filename = 'AvOwork/output/deltaGprediction/mousezfcellsCorr_qual_heatmap_withAnno_final_allzeros.pdf', height = 12, width = 14)


rowHist <- function(row, ptOrdGcTempl){
  y=as.numeric(ptOrdGcTempl[row,-1])
  par(mfrow=c(1,2), mar=c(5, 4, 4, 2) - 0.1)
  hist(y, xlab='predictedDg', breaks=15, main = ptOrdGcTempl[row,1])
  abline(v = mean(y, na.rm = T), col='red', xpd=F)
  legend('top', 'mean', lty=1, col = 'red', bty='n', xpd=T, inset = c(0,-0.05))
  qqnorm(y)
  qqline(y, xpd=F)
  test = shapiro.test(y)
  legend('topleft', paste0('SW test p-val = ', round(test$p.value, digits = 4)), cex = 0.7, bty='n')
}

diag <- take.diagonal(logPtPairs)
sapply(seq(1:10), function(x) rowHist(x,diag))

avgDeltaG <- function(ptPairs, no.reads){
  ptWithDg <- ptPairs %>%
    select(c('TemplatePrimer_pair', no.reads$cell[no.reads$num.reads > 40000])) %>%
    mutate(predictedDg = apply(.[,-1],1, mean, na.rm=TRUE),
           sd = apply(.[,-1],1, sd, na.rm=TRUE),
           n = is.na(.[,-1])) # Number of non NA va
  predictedDg <- ptWithDg[,c('TemplatePrimer_pair', 'predictedDg', 'sd', 'n')]
  return(predictedDg)
}
 
### Is the goodness of fit relative to the T?
load.cellAbundance <- function(cellAbFile) {
  sampleName <- gsub(gsub(cellAbFile, pattern = '.cellAbundance.+', replacement = ""), pattern = ".+/", replacement = '')
  cellAbundance <- read.csv(cellAbFile, header=TRUE, row.names=1)
  colnames(cellAbundance) <- c(paste0(sampleName,'.cell', gsub(colnames(cellAbundance), pattern = 'X', replacement = '')))
  return(cellAbundance)  
}


cellAbFile="~/mnt/edann/hexamers/rnaseq/mouse/testing/SvdB11d1-MitoTrackerThird-Satellites-Adult.cellAbundance.noN.csv"
cellAbundance<- load.cellAbundance(cellAbFile)

cellAb.cells <- cellAbundance[,colnames(corr)]
transcr.corr <- cor(cellAb.cells)
cors <- data.frame(PCC_Abundance=as.numeric(transcr.corr), PCC_DeltaG=as.numeric(corr))
ggplot(cors,aes(PCC_Abundance, PCC_DeltaG) ) + 
  geom_hex(binwidth=c(0.02,0.02)) +
  theme_classic()

ggsave("~/AvOwork/output/deltaGprediction/PCC_abundanceVSDeltaG.pdf")

