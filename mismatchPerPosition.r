## MISMATCH PER READ ANALYSIS
library(dplyr)
library(data.table)
library(ggplot2)

tab <- read.csv('~/mnt/edann/hexamers/mismatch/L1_mismatchfreq.csv', header = 1, row.names = 1)
tab <- tab[!grepl('N', rownames(tab)),]
melt(t(tab)) %>% ggplot(., aes(Var1, value, fill=Var2)) + geom_bar(stat='identity')

freqTab <- apply(tab, 2, function(x) x/sum(x, na.rm = TRUE))
melt(t(freqTab)) %>% ggplot(., aes(Var1, value, fill=Var2)) + geom_bar(stat='identity')
