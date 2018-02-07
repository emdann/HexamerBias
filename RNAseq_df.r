library(dplyr)
library(data.table)
library(ggplot2)
library(seqLogo)

dfHexMouserna <- fread('~/mnt/edann/hexamers/rnaseq/SvdB11d_sum.hexTab.wAb.csv', sep=',')
dfHexZFrna <- fread('~/mnt/edann/hexamers/rnaseq/gk2a-2.hexTable.wAb.csv', sep=',')


dfHexZFrna %>% mutate(label = ifelse(bindingEventsTempl>15000 | bindingEventsPrimer > 15000, V1, '')) %>%
  ggplot(., aes(mmEventsPrimer/bindingEventsPrimer, avgMmPrimer,color=sapply(V1, GCcont), label=V1)) + 
  geom_point() + #geom_text(nudge_y = 0.001, cex=2) +
  scale_color_gradient2(midpoint = 0.5) +
  labs(color='GC content') 
  annotate("text", label = paste0("PCC = ", round(cor(dfHexZFrna$bindingEventsTempl, dfHexZFrna$bindingEventsPrimer), digits = 2)), x = 5000, y = 50000, size = 4, colour = "black") +
  ggsave('~/AvOwork/output/RNAmodel/bindingEv_cor_GCcont.pdf')

pdf('AvOwork/output/SvdB11d4_pairs.pdf', height = 10, width = 10)
pairs(dfHexZFrna[,-1])
dev.off()


cor(dfHexZFrna$avgMmPrimer, log(dfHexZFrna$bindingEventsPrimer))

lm.fit = lm(formula = bindingEventsTempl ~ ., data=dfHexZFrna[,-1])

dfHexZFrna %>% mutate(binding_abundance_ratio=bindingEventsTempl/abundance) %>%
  mutate(label = ifelse(bindingEventsTempl>15000, V1, '')) %>% arrange(desc(bindingEventsTempl))
  ggplot(., aes( avgMmPrimer,bindingEventsPrimer, color=binding_abundance_ratio, label=label)) + geom_point() + 
  geom_text(cex=3, nudge_y = 0.1) + 
  scale_color_gradient2(midpoint = 1) +
  ggsave('~/AvOwork/output/RNAmodel/bigBinders_lowMM_primer.pdf')

files = list.files('mnt/edann/hexamers/rnaseq', pattern = '*primerDanRer.pwm.csv', full.names = TRUE)
seqs = gsub(gsub(files, pattern = '.+/', replacement = ''), pattern = '\\..+', replacement = '')
pwms = lapply(files, read.csv, row.names = 1)
names(pwms) <- seqs
ord_seq <- dfHexZFrna %>% filter(V1 %in% seqs) %>% mutate(GCcont = sapply(V1, GCcont)) %>% arrange(bindingEventsPrimer) %>% .$V1

pdf('AvOwork/output/RNAmodel/seqLogos_primer_ordBindingEvPrimer.pdf', onefile=T, height=10, width = 7)
ncol <- 1
nrow <- 4
s <- 1
while (s <= length(pwms)) {
  grid.newpage()
    for(row in 1:nrow){
    for(col in 1:ncol){
      vp <- viewport(x = (col-1)/ncol, y = 1-(row-1)/nrow, w
                     = 1/ncol, h =
                       1/nrow, just = c("left", "top"))
      pushViewport(vp)
      pwm = pwms[[ord_seq[s]]]
      mySeqLogo(pwm)
      grid.text(sprintf(ord_seq[s], row, col),
                x=0.5, y=0.9,
                just="top")
      upViewport()
      s <- s+1
    }
  }
}
dev.off()