### BS DEGRADATION IMPACT AND STRAND BIAS
library(dplyr)
library(data.table)
library(ggplot2)

cov.nuc <- fread("~/mnt/edann/BS_degradation/sorted_L1_trim1_R1_bismark_bt2_pe.highcov.windows.ratio.bed")
colnames(cov.nuc)[1:7] <- c('chr', 'start', 'end', 'id', 'coverage', 'strand', 'plus_minus_ratio')
colnames(cov.nuc) <- gsub(colnames(cov.nuc), pattern="[1234567890]+_", replacement = '')

# Is there a strand bias for most covered regions?
ggplot(cov.nuc, aes(plus_minus_ratio)) + geom_histogram(bins = 50) +
  theme_classic() +
  labs(x='log(+/-)', title='Strand ratio VAN1667 L1 (Most covered regions, sliding windows of 100 bps)') +
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size=14))


ggsave("~/AvOwork/output/strand_bias/strandratio_highcov_sortedL1_VAN1667.pdf")

# Are C rich regions subject to strand bias?
ggplot(cov.nuc, aes(x=num_C, y=plus_minus_ratio)) +
  # facet_wrap(~chr) +
  geom_point(alpha=0.3) +
  theme_classic() +
  xlab("% C") + ylab('log(+/-)') +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size=14))
    
ggsave("~/AvOwork/output/strand_bias/strand_biac_percC_L1_VAN1667.pdf", height = 10, width = 10)

