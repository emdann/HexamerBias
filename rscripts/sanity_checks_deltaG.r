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
    tab <- fread(paste('zcat',file), sep = ',', header=TRUE)
  }else{
    tab <- fread(file, sep = ',', header=TRUE) 
  }
p <- as.matrix(tab[,-1])
rownames(p)<- tab$V1
logMat <- log(p)
return(logMat)
  }

getDiag <- function(predDg){
  diag <- data.frame(predicted.deltaG = sapply(rownames(predDg), function(hex) predDg[hex,hex]))
  return(diag)
  }


## Concordance of NN model and diagonal
predDg <- read.csv(gzfile("~/mnt/edann/hexamers/VAN1667prediction/sorted_L1_trim1_R1_bismark_bt2_pe_ptDg_qual.csv.gz"), header = TRUE, row.names = 1)
diag <- getDiag(predDg)
## Load tabulated free energy
tabDg <- read.delim(gzfile("~/mnt/edann/hexamers/rand_hex_deltaG_ions.txt.gz"), sep=' ', header = FALSE, row.names = 1)
colnames(tabDg) <- 'freeEn'

dg.df <- cbind(pred=diag[match(rownames(tabDg), rownames(diag)),], tabDg)
dg.df <- dg.df %>% mutate(pred=log(pred), template=rownames(dg.df)) 
  # filter(pred==-Inf) %>%
dg.df %>%
  filter(!grepl('C|G', template)) %>% 
  ggplot(.,aes(pred,freeEn, label=template)) + geom_point() +
  # geom_smooth(method = lm, se = FALSE) +
  xlab('Predicted DeltaG') + ylab('NN DeltaG') +
  theme(axis.title = element_text(size = 20)) +
  # geom_text_repel() +
  xlim(-15,-3) +
  ylim(min(dg.df$freeEn), max(dg.df$freeEn))

## Simmetry
noCG.dg <- predDg[!grepl('C|G', rownames(predDg)), !grepl('C|G', colnames(predDg))]
noCG.dg <- noCG.dg[,match(row.names(noCG.dg), colnames(noCG.dg))] # Same order in cols and rows
noCg.dg.log <- log(noCG.dg)
noCg.dg.log[noCg.dg.log==-Inf]<-0
long.dg.log <- noCg.dg.log %>% mutate(template=rownames(noCg.dg.log)) %>% melt(value.name = 'predicted.Dg', variable.name = 'primer') 
long.dg.log$primer <- factor(long.dg.log$primer, levels = levels(factor(long.dg.log$template)) )
ggplot(long.dg.log, aes(template,primer,fill=predicted.Dg)) + geom_tile() +
  scale_fill_gradient2(midpoint = -10) +
  # scale_fill_gradientn(colours = brewer.pal(n = 9,'RdYlBu') +
  theme(axis.text.x = element_text(angle=90), axis.title = element_text(size = 20))

## Here I was trying to see clusters
pheatmap(noCg.dg.log, cluster_rows = TRUE, cluster_cols = FALSE, breaks = breaks, color = colorRampPalette(brewer.pal(n = 7, name =
                                                                                                                        "YlGnBu"))(100))
