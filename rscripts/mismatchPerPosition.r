## MISMATCH PER READ ANALYSIS
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)

## Compare primed region and used hexamer
counts=read.csv('mnt/edann/hexamers/rnaseq/hexUsage.csv')


counts %>% filter(primed_region_count < 20000 & usage_count < 20000) %>% ggplot(., aes(primed_region_count, usage_count, label=seq)) +
  geom_text(cex=2) +
  xlab('Usage as primed region') +
  ylab('Usage as hexamer') +
  ggsave('~/AvOwork/output/hexVSprimed_usage_low.pdf')

avgMM=read.csv('mnt/edann/hexamers/rnaseq/hexMMusage.csv')

df <- counts %>% mutate(avgMM=avgMM[match(counts$X, avgMM$X),]$avgMM) %>%
  filter(!grepl('N', X)) %>% mutate(gc=sapply(df$X, GCcont))
  %>% mutate(gc=ifelse(gc>=0.5, 'high', 'low'))

ggplot(df, aes(avgMM, usage, color=gc, label=X)) +
  geom_text(cex=2) +
  ylab('Usage as primer') +
  xlab('avg # mismatches') + 
  scale_color_gradient2(midpoint = 0.5) +
  labs(color='GC content') +
  ggsave('~/AvOwork/output/avgMMVShex_usage_gccont_rna.pdf')
  
ggplot(df, aes(x=avgMM, fill=gc)) +
  geom_histogram(binwidth=.2, alpha=.5, position="identity") +
  labs(fill="GC cont") +
  ggsave('~/AvOwork/output/avgMMdistGC_rna.pdf')

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

ggplot(d, aes(pos, freq)) + 
  geom_bar(stat='identity',aes(fill=MMtype), size=0, width = 1) +
  labs(x='Position', y='# reads') +
  geom_label(x=15, y=150000, label='Tot. no. of reads = 8245528', cex=5) +
  theme(axis.title = element_text(size = 20), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  ggsave(filename = '~/AvOwork/output/mismatch_per_position_abs_L1.pdf')

# n= 4003987
  
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
  