## MISMATCH PER READ ANALYSIS
library(dplyr)
library(data.table)
library(ggplot2)

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

  
ggplot(d, aes(pos, freq, fill=MMtype)) + 
  geom_bar(stat='identity', size=0, width = 1.2) +
  labs(x='Position') +
  ggsave(filename = '~/AvOwork/output/mismatch_per_position_abs_L1.pdf')
