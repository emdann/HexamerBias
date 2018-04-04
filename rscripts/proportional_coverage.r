### Proportional coverage ###
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)

prop.cov <- read.csv('~/mnt/edann/hexamers/VAN1667prediction/proportional_cov.csv', row.names = 1)
abundance <- read.csv(gzfile('~/mnt/edann/hexamers/VAN1667prediction/mm10.cellAbundance.noN.csv.gz'), row.names=1, header=FALSE)

df <- data.frame(prop.cov, abundance=abundance$V2[match(rownames(prop.cov), rownames(abundance))])


pred <- read.delim('~/mnt/edann/hexamers/VAN1667prediction/predictedCoverage_avgVAN1667.txt', row.names = 1)
df2 <- data.frame(prop.cov, pred.cov=pred$exp[match(rownames(prop.cov), rownames(pred))])
df2 %>% mutate(template= rownames(df2)) %>% ggplot(., aes(pred.cov, exp, label=template)) + geom_text()
