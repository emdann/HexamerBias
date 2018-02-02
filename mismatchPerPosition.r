## MISMATCH PER READ ANALYSIS
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)

## Compare primed region and used hexamer
counts=read.delim('mnt/edann/hexamers/mismatch/L1_R1_hexVSprimed.count.txt')

counts %>% filter(primed_region_count < 20000 & usage_count < 20000) %>% ggplot(., aes(primed_region_count, usage_count, label=seq)) +
  geom_text(cex=2) +
  xlab('Usage as primed region') +
  ylab('Usage as hexamer') +
  ggsave('~/AvOwork/output/hexVSprimed_usage_low.pdf')

avgMMtol=read.delim('mnt/edann/hexamers/mismatch/L1_R1_avgMM_tolerateBSmm.txt', header=FALSE)
  
counts %>% mutate(avgMM=avgMMtol[match(counts$seq, avgMMtol$V1),]$V2) %>%
  ggplot(., aes(avgMM, primed_region_count, label=seq)) +
  geom_text(cex=2) +
  ylab('Usage as primed region') +
  xlab('avg # mismatches') +
  ggsave('~/AvOwork/output/avgMMVShex_usage_tolBS.pdf')
  
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
mmmatrix <- fread("~/mnt/edann/hexamers/mismatch/L1_R1_primed_seq.original.csv", header = TRUE)
colnames(mmmatrix)[1]<- 'hex'
normmatrix <- apply(mmmatrix[,-1], 2, function(x) x/sum(x))
rownames(normmatrix)<-mmmatrix$primedSeq
mat_hm<- function(ixs, cluster_rows = FALSE, cluster_cols = FALSE){
  # drows = cor(normmatrix[ixs,ixs], method = 'euclidean')
  # dcols = dist(t(normmatrix[ixs,ixs]), method = 'euclidean')
  pheatmap(normmatrix[ixs,ixs], cluster_rows = cluster_rows, clustering_distance_rows = 'correlation' , cluster_cols = cluster_cols , clustering_distance_cols = 'correlation', border_color = NA, labels_row = mmmatrix[ixs,]$primedSeq)
}
drows = dist(normmatrix, method = 'euclidean')
dcols = dist(t(test), method = "euclidean")
  