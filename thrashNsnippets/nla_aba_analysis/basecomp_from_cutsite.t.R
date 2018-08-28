library(reshape2)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)

file="~/mnt/edann/AbaNla/CG-abaBS-1000ES-3n1_AHFKWYBGX5_S13_L004_R1_001.basematrix.csv"
load.basematrix <- function(file) {
  base.mat <- read.csv(file)
  colnames(base.mat) <- c('base', gsub(colnames(base.mat), pattern = 'X', replacement = '')[-1])
  long.mat <- melt(base.mat, id.vars = 'base', variable.name = "position")
  long.mat <- filter(long.mat, long.mat$base!='N')
  return(long.mat)
}

tabs <- lapply(list.files("~/mnt/edann/AbaNla/", full.names = TRUE), load.basematrix)
t <- tabs[[1]]
n <- 1
for (tab in tabs[-1]) {
  t <- merge(t, tab, by=c('position', 'base'), suffixes = c(n, n+1))  
  n <- n+1
}

t <- t %>% mutate(sum = value1+value2+value3+value4)
p <- ggplot(filter(t, as.numeric(position) < 33), aes(position, sum, group=base, color=base)) + 
  theme_classic() +
  geom_line() 

p + xlab('Distance from cut site (bps)') + ylab("Frequency") +
  scale_color_discrete(name="") +
  ggsave("~/AvOwork/output/Aba_Nla_analysis/basecomp_from_aba_cutsite.pdf")
