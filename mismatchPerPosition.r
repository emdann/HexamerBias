## MISMATCH PER READ ANALYSIS
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)

## Mismatch per position
tab <- read.csv('~/mnt/edann/hexamers/mismatch/L1_mismatchfreq.csv', header = 1, row.names = 1)
tab <- tab[!grepl('N', rownames(tab)),]
melt(t(tab)) %>% ggplot(., aes(Var1, value, fill=Var2)) + geom_bar(stat='identity')

freqTab <- apply(tab, 2, function(x) x/sum(x, na.rm = TRUE))
d <- melt(t(tab)) %>% 
  rename(pos=Var1, MMtype=Var2, freq=value) %>%
  mutate(pos=as.numeric(gsub(pos, pattern = 'X', replacement = '')), freq=ifelse(is.na(freq),0,freq)) %>%
  arrange(pos)

d <- melt(t(freqTab)) %>% 
  rename(pos=Var1, MMtype=Var2, freq=value) %>%
  mutate(pos=as.numeric(gsub(pos, pattern = 'X', replacement = '')), freq=ifelse(is.na(freq),0,freq)) %>%
  arrange(pos) 

ggplot(d, aes(pos, freq, fill=MMtype)) + 
  geom_bar(stat='identity', size=0, width = 1.2) +
  labs(x='Position') +
  ggsave(filename = '~/AvOwork/output/mismatch_per_position_abs_L1.pdf')

## Mismatch matrix
mmmatrix <- fread("~/HexamerBias/test_mmMatrix.csv", header = TRUE)
colnames(mmmatrix)[1]<- 'primedSeq'
normmatrix <- apply(mmmatrix[,-1], 2, function(x) x/sum(x))
mat_hm<- function(ixs){
  # drows = cor(normmatrix[ixs,ixs], method = 'euclidean')
  # dcols = dist(t(normmatrix[ixs,ixs]), method = 'euclidean')
  pheatmap(normmatrix[ixs,ixs], cluster_rows = TRUE, clustering_distance_rows = 'correlation' , cluster_cols = TRUE , clustering_distance_cols = 'correlation', border_color = NA, labels_row = mmmatrix[ixs,]$primedSeq)
}
drows = dist(normmatrix, method = 'euclidean')
dcols = dist(t(test), method = "euclidean")
  