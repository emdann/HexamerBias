## BS degradation: methylation bias
library(dplyr)
library(ggplot2)
library(data.table)

# setwd('/hpc/hub_oudenaarden/edann/BS_degradation')


# highcov.met1 <- fread("/hpc/hub_oudenaarden/edann/crypts_bs/VAN1667/L1_trim1_R1_bismark_bt2_pe.deduplicated.bismark.thresh5.cov.bed")
# highcov.met2 <- fread("/hpc/hub_oudenaarden/edann/hexamers/kaester/met_extraction/ERR454965_1_val_1_bismark_bt2.deduplicated.bismark.thresh5.cov")

highcov.met1 <- fread("~/mnt/edann/crypts_bs/VAN1667/L1_trim1_R1_bismark_bt2_pe.deduplicated.bismark.thresh5.cov.bed")
highcov.met2 <- fread("~/mnt/edann/hexamers/kaester/met_extraction/ERR454965_1_val_1_bismark_bt2.deduplicated.bismark.thresh5.cov")

colors <- c('ampPBAT', 'BS-seq')

ggplot(data=NULL) +
  geom_density(data=highcov.met1, aes(V4, color=colors[1], fill=colors[1]), alpha=0.5) +
  geom_density(data=highcov.met2, aes(V4, color=colors[2], fill=colors[2]), alpha=0.5) +
  theme_classic() +
  xlab('Methylation fraction') +
  scale_color_discrete(name="Protocol", labels=colors) +
  scale_fill_discrete(name="Protocol", labels=colors) +
  ggtitle('Meth fraction mouse crypts (coverage > 5)') +
  theme(axis.title = element_text(size = 18), title = element_text(size=20)) +
  ggsave('~/AvOwork/output/meth_frac_pbatVSbsseq.pdf')
