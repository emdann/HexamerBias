## COVERAGE PEAKS 
library(data.table)
library(dplyr)
library(ggplot2)

peakAnn <- fread("mnt/edann/coverage_peaks/multipeaks.annotatePeaks.homer")
randAnn <- fread("mnt/edann/coverage_peaks/multipeaks.random.annotatePeaks.homer")

pdf("AvOwork/output/covPeaks_distTSS_boxplot.pdf")
boxplot(abs(peakAnn$`Distance to TSS`), abs(randAnn$`Distance to TSS`), outline = FALSE, varwidth = TRUE, names = c("Coverage peaks", "random"), ylab='Distance to TSS')
dev.off()

g <- peakAnn %>% mutate(Annotation, gsub(pattern = " \\(.+\\)",replacement = "", x=Annotation)) %>%
  ggplot(., aes(Annotation)) + geom_bar()
