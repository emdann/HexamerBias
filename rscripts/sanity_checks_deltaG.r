## Sanity checks on predicted DeltaG
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)


loadDgMatrix <- function(file, compression='gzip'){
  if (compression=='gzip') {
    # tab <- fread(paste('zcat',file), sep = ',', header=TRUE)
    tab <- read_csv(file) ## <---- WON'T KEEP ROWNAMES!
  }else{
    tab <- fread(file, sep = ',', header=TRUE) 
  }
  p <- as.matrix(tab[,-1])
  rownames(p)<- tab$V1
  logMat <- log(p)
  return(logMat)
  }

make_pair_df <- function(dgMat){
  pairDf <- dgMat %>% 
    melt(value.name = 'dG') %>% 
    mutate(ptPair=paste0(Var1,'.',Var2)) %>% 
    select(ptPair, dG)
  return(pairDf)
}


getDiag <- function(predDg){
  diag <- data.frame(predicted.deltaG = sapply(rownames(predDg), function(hex) predDg[hex,hex]))
  return(diag)
  }

getDiag.pair <- function(Dg.pair){
  diag <- Dg.pair[substr(Dg.pair$ptPair,1,6)==substr(Dg.pair$ptPair,8,13),]
  return(diag)
}

# ## Concordance of NN model and diagonal
# predDg <- read.csv("~/mnt/edann/hexamers/strand_specific/VAN1667_se_ptDg_qual.csv", header = TRUE, row.names = 1)
# diag <- getDiag(predDg)
# ## Load tabulated free energy
# tabDg <- read.delim(gzfile("~/mnt/edann/hexamers/rand_hex_deltaG_ions.txt.gz"), sep=' ', header = FALSE, row.names = 1)
# colnames(tabDg) <- 'freeEn'
# 
# dg.df <- cbind(pred=diag[match(rownames(tabDg), rownames(diag)),], tabDg)
# dg.df <- dg.df %>% mutate(pred=log(pred), template=rownames(dg.df)) 
#   # filter(pred==-Inf) %>%
# dg.df.clean <- dg.df %>% 
#   # filter(!grepl('C', template)) %>% 
#   filter(pred!='-Inf') %>%
#   mutate(lab=ifelse(freeEn< -5 & -pred <10 | freeEn> -5 & -pred >10.5, NA, template))
# 
# ggplot(dg.df.clean,aes(freeEn, -pred, labels=lab)) + 
#   geom_point(size=3) + 
#   # stat_dens2d_filter(geom='text_repel', keep.fraction=0.1) +
#   ylab('Predicted DeltaG') + xlab('NN DeltaG') +
#   theme_bw() 
#   theme(axis.title = element_text(size = 30), 
#         axis.text = element_text(size=25), 
#         title = element_text(size=22)) 
#   # ylim(min(dg.df$freeEn), max(dg.df$pred)) +
#   xlim(min(dg.df$freeEn), max(dg.df$freeEn))
#   # ggtitle('Strand specific DeltaG prediction')
# 
# ggsave('~/AvOwork/output/deltaGprediction/sanity_checks/predBsVSnnmodel_noC_strandSpecific.pdf')
# 
# ## Simmetry
# noCG.dg <- predDg[!grepl('C|G', rownames(predDg)), !grepl('C|G', colnames(predDg))]
# noCG.dg <- noCG.dg[,match(row.names(noCG.dg), colnames(noCG.dg))] # Same order in cols and rows
# noCg.dg.log <- log(noCG.dg)
# noCg.dg.log[noCg.dg.log==-Inf]<- -18
# long.dg.log <- noCg.dg.log %>% mutate(template=rownames(noCg.dg.log)) %>% melt(value.name = 'predicted.Dg', variable.name = 'primer') 
# long.dg.log$primer <- factor(long.dg.log$primer, levels = levels(factor(long.dg.log$template)) )
# ggplot(long.dg.log, aes(template,primer,fill=predicted.Dg)) + geom_tile() +
#   scale_fill_gradient2(midpoint = -10) +
#   # scale_fill_gradientn(colours = brewer.pal(n = 9,'RdYlBu') +
#   theme(axis.text.x = element_text(angle=90), axis.title = element_text(size = 20))
# 
# ## Here I was trying to see clusters
# noCg.dg.log.srt <- noCg.dg.log[sort(rownames(noCg.dg.log)),]
# pheatmap(noCg.dg.log, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(brewer.pal(n = 7, name =
#                                                                                                                         "YlGnBu"))(100))
# ## Assuming conversion
# CH.hex <- rownames(predDg)[grepl(rownames(predDg), pattern='C[^G]')]
# bs.convert <- function(seq){  return(gsub(seq, pattern = 'C', replacement = 'T')) }
# 
# df <- data.frame(CH.hex, non.conv.Dg = sapply(CH.hex, function(seq) predDg[seq, seq]),conv.Dg = sapply(CH.hex, function(seq) predDg[seq, bs.convert(seq)]), nnDg=sapply(CH.hex, function(seq) tabDg[bs.convert(seq),]))
# ggplot(df, aes(log(conv.Dg), nnDg)) + geom_point() +
#   xlab('deltaG CH template non conv') +
#   ylab('NN delta G conv primer') +
#   theme(axis.title = element_text(size = 18))
# 
# 
# ggsave("~/AvOwork/output/deltaGprediction/sanity_checks/non_conversion_CH.pdf")




