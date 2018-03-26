library(dplyr)
library(ggplot2)

minus <- read.delim("~/mnt/edann/hexamers/VAN1667prediction/mm10.minusstrand.abundance.noN.txt", sep='\t')
plus <- read.csv(gzfile("~/mnt/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz"), header = F, col.names = c('abundance','kmer'))

d <- merge(minus, plus, by='abundance')
ggplot(d, aes(kmer.x, kmer.y)) + geom_point() +
  xlab('T abundance reverse strand') +
  ylab('T abundance forward strand') +
  theme(axis.title = element_text(size = 18))
