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

freqTab <- apply(tab, 2, function(x) x/8245528)
d <- melt(t(tab)) %>% 
  rename(pos=Var1, MMtype=Var2, freq=value) %>%
  mutate(pos=as.numeric(gsub(pos, pattern = 'X', replacement = '')), freq=ifelse(is.na(freq),0,freq)) %>%
  arrange(pos)

d <- melt(t(freqTab)) %>% 
  rename(pos=Var1, MMtype=Var2, freq=value) %>%
  mutate(pos=as.numeric(gsub(pos, pattern = 'X', replacement = '')), freq=ifelse(is.na(freq),0,freq)) %>%
  arrange(pos) 

d1 <- filter(d, !grepl(d$MMtype, pattern = 'N->.+')) %>%
  filter(!MMtype %in% c('A->A', 'G->G', 'C->C', 'T->T')) %>%
  filter(pos< 60) 
  
p <- ggplot(d1, aes(pos, freq)) + 
    theme_classic() +
    geom_bar(stat='identity', size=0, width = 0.7) +
    labs(x='Position', y='# reads') 
p + 
  geom_label(x=40, y= 120000, label='Tot. no. of reads = 8245528', cex=10) +
    # scale_fill_discrete(name='Mismatch') +
    theme(axis.title = element_text(size = 30), 
          axis.text = element_text(size=25), 
          title = element_text(size=22),
          legend.title = element_blank(),
          legend.text = element_text(size=25), 
          legend.key.size = unit(1,"cm"))

ggsave(filename = '~/AvOwork/formatted_figs/Mismatch figure//mismatch_per_position_abs_L1.pdf')

templ.base <- d1 %>% 
  filter(pos<7) %>%
  mutate(original.base=substr(MMtype, start = 1, stop = 1),
         obs.base = substr(MMtype, start = 4, stop = 4)) %>%
  ggplot(., aes(pos, freq*100, fill=original.base)) + 
  theme_classic() +
  scale_x_continuous(breaks=seq(1,6)) +
  geom_bar(stat='identity', position='fill',size=0, width = 0.7, ) +
  labs(x='Position', y='frac. of reads') +
  theme(axis.title = element_text(size = 30), 
        axis.text = element_text(size=25), 
        title = element_text(size=22),
        legend.title = element_blank(),
        legend.text = element_text(size=25), 
        legend.key.size = unit(1,"cm")
        )
  NULL

primer.base <- d1 %>% 
  filter(pos<7) %>%
  mutate(original.base=sapply(substr(MMtype, start = 1, stop = 1), rev.comp),
         obs.base = substr(MMtype, start = 4, stop = 4)) %>%
  ggplot(., aes(pos, freq, fill=obs.base)) + 
  theme_classic() +
  scale_x_continuous(breaks=seq(1,6)) +
  geom_bar(stat='identity', position='fill',size=0, width = 0.7, ) +
  labs(x='Position', y='frac. of reads') +
  theme(axis.title = element_text(size = 30), 
        axis.text = element_text(size=25), 
        title = element_text(size=22),
        legend.title = element_blank(),
        legend.text = element_text(size=25), 
        legend.key.size = unit(1,"cm")
  )
  NULL

primer.base + ggsave('~/AvOwork/formatted_figs/Mismatch figure/primer_mm.pdf')
templ.base + ggsave('~/AvOwork/formatted_figs/Mismatch figure/templ_mm.pdf')
# n= 4003987

d1 %>% 
  filter(pos<7) %>%
  mutate(orig.base=substr(MMtype, start = 1, stop = 1),
       obs.base = substr(MMtype, start = 4, stop = 4)) %>%
  ggplot(., aes(y=freq*100, fill=MMtype)) + 
  theme_bw() +
  geom_bar(aes(x=pos),stat='identity', position='dodge',size=0) +
  facet_grid(orig.base~obs.base, scales='free_x', labeller=label_both) +
  scale_x_continuous(breaks=seq(1,6)) +
  scale_fill_brewer(palette='Set3') +
  ylab('% reads') + xlab("Position") +
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size=20), 
        title = element_text(size=22),
        strip.text = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=25), 
        legend.key.size = unit(1,"cm")
        ) 
  ggsave("~/AvOwork/formatted_figs/Mismatch figure/primer_template_mm.pdf")
 
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
  