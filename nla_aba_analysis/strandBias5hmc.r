library(data.table)
library(dplyr)
library(plyr)

bot_5hmc <- fread("~/mnt/edann/AbaNla/CG-abaBS-1000ES-3n1/CG-abaBS-1000ES-3n1_lmerged_R1.CTOB.uniq5hmc.bed", col.names = c('chr', 'start', 'end'))
bot <- bot_5hmc %>% mutate(strand='bot')
top_5hmc <- fread("~/mnt/edann/AbaNla/CG-abaBS-1000ES-3n1/CG-abaBS-1000ES-3n1_lmerged_R1.CTOT.uniq5hmc.bed", col.names = c('chr', 'start', 'end'))
top <- top_5hmc %>% mutate(strand='top')

hmc <- rbind(bot,top)

plotStrandHmc <- function(chrom){
  d <- hmc %>% filter(chr==chrom)
  ggplot(d, aes(x=start, fill=strand)) + 
    geom_histogram(data = subset(d, strand %in% 'top'), aes(x = start, y = ..density..),  binwidth = 50000) +
    geom_histogram(data = subset(d, strand %in% 'bot'), aes(y = -..density..) , binwidth = 50000) +
    xlab(chrom)
}

pdf("~/AvOwork/output/Aba_Nla_analysis/strand_distribution_5hmc_density.pdf", onefile = TRUE, width = 12)
for (chr in unique(hmc$chr)) {
  p <- plotStrandHmc(chr)
  plot(p)
}
dev.off()
