### TEST COVERAGE PREDICTION ACCURACY ###
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
options(stringsAsFactors = FALSE)

## COMPUTE RSE
res.std.err <- function(obs,exp){ # Make sure templates are matched
  rse <- sqrt((1/length(obs))*sum((obs-exp)^2))
  return(rse)
}

pred <- pred[match(rownames(obs), rownames(pred)),]
pred.dataframe <- obs %>% mutate(exp=pred$exp, err=pred$err)
rse <- sapply(obs, function(x) res.std.err(x, pred$exp))
nreads <- colSums(obs)
df <- data.frame(sample = names(rse), n.reads=abs(nreads-mean(nreads[-1])), rse=rse)
ggplot(data = df, aes(y=rse,x=n.reads, label=sample)) + 
  theme_classic() +
  geom_point() + geom_text_repel(cex=3) +
  xlab("No. of reads - avg no. reads training set") +
  ylab("Residual Standard Error") +
  theme(axis.title = element_text(size = 18))
ggsave("~/AvOwork/output/deltaGprediction/deltaNumberOfReads.pdf")

## RSE for different predicted coverages
pred <- read.csv("~/mnt/edann/hexamers/zf_prediction/zf_covPrediction_avgOnNreads.csv", row.names = 1)
obs <- obs.zf[,grepl(colnames(obs.zf), pattern='R1')]

plot.rse.for.predictions <- function(pred.tab, obs.tab){
  ## Are the rownames matched?
  pred <- pred[match(rownames(obs), rownames(pred)),]
  rse.nreads <- sapply(obs, function(sample) sapply(pred, function(predDg) res.std.err(sample, predDg)))
  long.rse.nreads <- melt(rse.nreads, varnames = c('training.set', 'sample'), value.name = 'RSE') %>%
    mutate(sample=gsub(sample, pattern='_trim2|_bismark.+', replacement = ''))
    # mutate(num.reads = num.reads[sample], training.set=as.numeric(substr(training.set,4,100)))
    
    ggplot(long.rse.nreads, aes(training.set, RSE, group=sample, color=)) + 
    theme_classic() +
    geom_point() + geom_line() +
    ylab('Residual Standard Error') + xlab("# of averaged samples") +
    theme(axis.title = element_text(size = 20), legend.title = element_text(size=14), legend.text = element_text(size=10))
  
}

ggsave("~/AvOwork/output/deltaGprediction/no.reads_in_avg_zf.pdf")

