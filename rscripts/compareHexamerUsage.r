### COMPARE HEXAMER USAGE ###
# Load from ptCounts tables
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

loadPtMatrix <- function(file, compression='gzip'){
  if (compression=='gzip') {
    tab <- fread(paste('zcat',file), sep = ',', header=TRUE)
  }else{
    tab <- fread(file, sep = ',', header=TRUE) 
  }
  p <- as.matrix(tab[,-1])
  rownames(p)<- tab$V1
  return(p)
}

pt.VAN1667 <- loadPtMatrix("~/mnt/edann/hexamers/strand_specific/VAN1667_se.ptCounts.qualFilt.parallel.csv", compression = 'no')
save(pt.VAN1667, file='~/AvOwork/ptCounts.strandSpec.VAN1667.RData')
pt.VAN2408 <- loadPtMatrix('~/mnt/edann/crypts_bs/VAN2408/CM1_tr2_R1_bismark_bt2.ptCounts.qualFilt.parallel.csv', compression='no')

## Compare hexamer profiles
make.hex.usage.df <- function(ptTab, type = 'primer', scale =TRUE){
  if (type=='primer') { hex.usage <- colSums(ptTab)  }
  if (type=='template') { hex.usage <- rowSums(ptTab)  }
  sample.name <- gsub(deparse(substitute(ptTab)), pattern = 'pt.', replacement = '')
  if (scale) { 
    hex.df <- data.frame(names(hex.usage), hex.usage/sum(hex.usage), row.names = NULL) 
  }else{
    hex.df <- data.frame(names(hex.usage), hex.usage, row.names = NULL)
    }
  colnames(hex.df) <- c(type,sample.name)
  return(hex.df)
}

VAN1667.primer.usage <- make.hex.usage.df(pt.VAN1667, type = 'primer', scale=TRUE)
VAN2408.primer.usage <- make.hex.usage.df(pt.VAN2408, type = 'primer', scale=TRUE)

primer.usage.df <- merge(VAN1667.primer.usage, VAN2408.primer.usage, by = 'primer')
d <- melt(primer.usage.df, variable.name = 'sample', value.name = 'usage')
ggplot(d, aes(x = primer, fill=sample)) + 
    theme_classic() +
    geom_col(data = filter(d, sample=='VAN1667'), aes(y = usage), position=position_dodge(0.9), width = 0.5) +
    geom_col(data = filter(d, sample=='VAN2408'), aes(y = -usage), alpha=1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank())

primer.usage.df %>% mutate(label = ifelse(VAN1667 > 0.005 | VAN2408 > 0.005, as.character(primer),'')) %>%
ggplot(., aes(VAN1667, VAN2408, label=label)) + 
  theme_classic() +
  geom_point() +
  geom_text_repel() +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14)) +
  xlab('Machine-mixed primers relative abundance') + ylab('Hand-mixed primers relative abundance')

ggsave('~/AvOwork/hexamer_usage_primer_hadmixed_scatter.pdf')

## Cumulative distribution of template sequences 

VAN1667.template.usage <- make.hex.usage.df(pt.VAN1667, type = 'template', scale=TRUE)
VAN2408.template.usage <- make.hex.usage.df(pt.VAN2408, type = 'template', scale=TRUE)
template.usage.df <- merge(VAN1667.template.usage, VAN2408.template.usage, by = 'template')
template.usage.df 
d <- data.frame(VAN1667.cumsum = cumsum(sort(template.usage.df$VAN1667, decreasing = TRUE)), 
                VAN2408.cumsum = cumsum(sort(template.usage.df$VAN2408, decreasing = TRUE)))
d1 <- d %>% mutate(ix=as.numeric(row.names(d))) %>%
  melt(variable.name='sample', id.vars=c('ix')) 
ggplot(d1, aes(ix,value, color=sample)) + geom_point() +   
  ggtitle('Cum distribution template sequences') +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14),
        title = element_text(size=22)) 
  
ggsave('~/AvOwork/cumsum_templates.pdf')

## Compare experiments from the same batch
MD <- loadPtMatrix("~/mnt/edann/crypts_bs/VAN2408/MD1_tr2_R1_bismark_bt2.ptCounts.qualFilt.parallel.csv", compression='no')
pt.MD <- MD

MD.primer.usage <- make.hex.usage.df(pt.MD, type='primer', scale=TRUE)

handMixed.df <- merge(primer.usage.df, MD.primer.usage) %>%
  mutate(CM = VAN2408) 

handMixed.df %>% mutate(label = ifelse(CM > 0.005 | MD > 0.005, as.character(primer),'')) %>%
  ggplot(., aes(CM, MD, label=label)) + 
  theme_classic() +
  geom_point() +
  geom_text_repel() +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14)) +
  xlab('CM1') + ylab('MD1')

ggsave('~/AvOwork/handMixed_primer_usage.pdf')

## My BS data
cele.noBS <- loadPtMatrix("~/mnt/edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-CeleTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv", compression = 'no')
zf.noBS <- loadPtMatrix("~/mnt/edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-zfishTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv", compression = 'no')

cele.primer.usage <- make.hex.usage.df(cele.noBS, type="primer", scale=TRUE)
zf.primer.usage <- make.hex.usage.df(zf.noBS, type="primer", scale=TRUE)

cele.template.usage <- make.hex.usage.df(cele.noBS, type="template", scale=TRUE)
cele.abs.template.usage <- make.hex.usage.df(cele.noBS, type="template", scale=F)
zf.template.usage <- make.hex.usage.df(zf.noBS, type="template", scale=T)
zf.abs.template.usage <- make.hex.usage.df(zf.noBS, type="template", scale=F)

merge(cele.template.usage, zf.template.usage, by='template') %>%
  mutate(label = ifelse(cele.noBS > 0.002 | zf.noBS > 0.003, as.character(template),'')) %>%
  ggplot(., aes(cele.noBS, zf.noBS, label=label)) + 
  theme_classic() +
  geom_point() +
  geom_text_repel() +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14)) +
  xlab('C.elegans primer usage') + ylab('D.rerio primer usage') +
  geom_abline(intercept=0, slope=1, color='red') + 
  ggsave("~/AvOwork/output/MyData/template_usage_celVSzf_noBS.pdf")


### Compare with total kmer abundance
cele.abundance <- read_csv("~/mnt/edann/hexamers/genomes_kmers/WBcel235.kmerAbundance.csv", col_names=FALSE)
colnames(cele.abundance) <- c('hex', 'abundance')


zf.abundance %>% 
  rename(template=hex) %>% 
  inner_join(., zf.template.usage, by='template') %>% 
  ggplot(., aes(abundance, zf.noBS)) + 
  geom_point(alpha=0.2) + 
  theme_classic() + 
  ylab('template usage') + xlab('template abundance') +
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size=20))
ggsave('~/AvOwork/output/MyData/abVStemplate_zf.pdf')

zf.abundance %>% 
  rename(template=hex) %>% 
  inner_join(., zf.abs.template.usage, by='template') %>% 
  ggplot(., aes(abundance, zf.noBS)) + 
  geom_point(alpha=0.2) + 
  theme_classic() + 
  ylab('template usage') + xlab('template abundance') +
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size=20)) +
  coord_cartesian(xlim=c(0,1000000), ylim=c(0,20000))

cele.abundance %>% 
  rename(template=hex) %>% 
  inner_join(., cele.abs.template.usage, by='template') %>% 
  ggplot(., aes(abundance, cele.noBS)) + 
  geom_point(alpha=0.2) + 
  theme_classic() + 
  ylab('template usage') + xlab('template abundance') +
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size=20)) +
  coord_cartesian(xlim=c(0,200000), ylim=c(0,20000))

## Try again original model
tabDg <- read.delim(gzfile("~/mnt/edann/hexamers/rand_hex_deltaG_ions.txt.gz"), sep=' ', header = FALSE, col.names=c('template', 'dG'))

zf.noBS.pair <- make_pair_df(zf.noBS)
zf.diag <- getDiag.pair(zf.noBS.pair)

zf.ab <- zf.abundance %>% 
  rename(template=hex) %>% 
  inner_join(., zf.abs.template.usage, by='template') %>%
  inner_join(., tabDg, by='template') 

zf.ab %>%
  filter(zf.noBS>500) %>%
  ggplot(., aes(dG, log(zf.noBS/abundance))) +
  geom_point(alpha=0.3)

zf.df <- zf.diag %>%
  mutate(hex=substr(ptPair,1,6)) %>%
  inner_join(., zf.abundance, by='hex') %>%
  rename(template=hex, cov=dG) %>%
  inner_join(., tabDg, by='template') 

zf.df %>%
  filter(cov>250) %>%
  ggplot(., aes(dG, log(cov/abundance))) +
  geom_point(alpha=0.3) +
  xlab('Delta G NN') + ylab('log(pt/T)') +
  ggtitle('Filtering pt < 250') +
  theme(axis.text = element_text(size=20),
      axis.title = element_text(size=24),
      title=element_text(size=28))

ggsave("~/AvOwork/output/MyData/model_fitting_zf_noBS_diag")


cors<-c()
for (n in seq(100,1900, by = 100)) {
  # print(n)
  pcc <- with(filter(zf.df, cov>n ),
       cor(dG,log(cov/abundance), use='pairwise.complete.obs'))
  cors <- c(cors, pcc)
}
plot(seq(100,1900,by=100), cors, type='b', ylim=c(0,-1))

cors.tot<-c()
for (n in seq(100,1900, by = 100)) {
  # print(n)
  pcc <- with(filter(zf.ab, zf.noBS>n ),
              cor(dG,log(zf.noBS/abundance), use='pairwise.complete.obs'))
  cors.tot <- c(cors.tot, pcc)
}

plot(seq(100,1900,by=100), cors.tot, type='b', ylim=c(0,-1))
