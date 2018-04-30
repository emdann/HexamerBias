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
BS <- loadPtMatrix("~/mnt/edann/noPreAmp_crypts/CG-pbat-1xPA-gDNA-crypts_lmerged_R1_val_1_bismark_bt2_pe.ptCounts.qualFilt.parallel.csv", compression = 'no')
pt.BS <- BS

BS.primer.usage <- make.hex.usage.df(pt.BS, type="primer", scale=TRUE)

merge(BS.primer.usage, MD.primer.usage, by='primer') %>%
  mutate(label = ifelse(BS > 0.005 | MD > 0.005, as.character(primer),'')) %>%
  ggplot(., aes(BS, MD, label=label)) + 
  theme_classic() +
  geom_point() +
  geom_text_repel() +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14)) +
  xlab('BS_Emma') + ylab('VAN2408_MD1')

ggsave("~/AvOwork/handMixed_primer_usage_meVSanna.pdf")

### Compare with total kmer abundance
abundance <- read.csv(gzfile('~/mnt/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz'), row.names=1, header=FALSE)
abundance <- abundance %>% mutate(template=rownames(abundance))
ab.primer.df <- base::merge(abundance,primer.usage.df, by='primer')
ab.template.df <- base::merge(abundance,template.usage.df, by='template')
