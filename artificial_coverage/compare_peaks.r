### Compare experimental and predicted coverage
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)
library(rtracklayer)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")
# biocLite("similaRpeak")
# biocLite('ChIPpeakAnno')
# biocLite('plyranges')
library(Gviz)
library(similaRpeak)
# library(plyranges)

# bed.intersect <- read.table('~/mnt/edann/hexamers/strand_specific/random.42.1000smp.chr1.intersect.bed')
# colnames(bed.intersect) <- c('chr', 'start', 'end', 'exp', 'pred')
# 
# bed.intersect %>% 
#   mutate(coverage=ifelse(exp >20, 20,exp)) %>%
#   mutate(coverage=cut(exp,c(0,0.0001,3,5,10,max(exp)), right=FALSE, include.lowest = TRUE)) %>%
#   ggplot(., aes(y=pred, x=coverage, group=coverage, fill=coverage)) +
#   # facet_wrap(~coverage) +
#   theme_classic() +
#   geom_boxplot(varwidth = TRUE, outlier.alpha=0.01, lwd=1.5) +
#   scale_fill_brewer() +
#   scale_x_discrete(labels=c('0', '< 3', '< 5', '< 10', '> 10')) +
#   xlab("experimental coverage") + ylab("Predicted coverage") +
#   theme(axis.title = element_text(size = 30), 
#         axis.text = element_text(size=25), 
#         title = element_text(size=22)) +
#   guides(fill='none')
  
## Trying to use rtracklayer

make.base.res.bw <- function(bw){
  base.res.bw.list <- tile(bw, width=1)
  base.res.bw <- base.res.bw.list@unlistData
  base.res.bw$score <- rep(exp.bw$score,width(exp.bw))
  return(base.res.bw)
}


make.predVSexp.range <- function(pred.bw, exp.bw){
  ranges <- subsetByOverlaps(exp.bw,pred.bw)
  hits <- findOverlaps(exp.bw, pred.bw)
  idx <- subjectHits(hits)
  values <- DataFrame(pred = pred.bw$score[idx])
  mcols(ranges) <- c(mcols(ranges), values)
  return(ranges)
} 

add.id <- function(subset.common){
  # print(subset.common@ranges@start)
  breaks <- sapply(seq_along(subset.common@ranges@start), function(ix) 
    subset.common@ranges@start[ix+1] - subset.common@ranges@start[ix]!=1)
  # print(breaks)
  breaks <- unlist(breaks)
  break.points <- unique(subset.common[which(breaks)+1]@ranges@start)
  ids <- cut(subset.common@ranges@start, breaks = c(0,break.points, max(break.points)+10000), include.lowest = TRUE)
  subset.common$id <- paste(subset.common@seqnames,ids)
  return(subset.common)
  }

## Subset
subsetByRegion <- function(bw, chrom, start, end){
  subset.bw <- bw[bw@seqnames==chrom &
                  bw@ranges@start >= start &
                  bw@ranges@start < end]
  return(subset.bw)  
}

### Normalizing by avg
normalize.coverage <- function(bw){
  bw.norm <- bw
  for (col in colnames(values(bw)[sapply(values(bw), is.numeric)])) {
    # bw.norm@elementMetadata[col][[1]] <- bw.norm@elementMetadata[col][[1]]/mean(bw.norm@elementMetadata[col][[1]])
    bw.norm@elementMetadata[col][[1]] <- bw.norm@elementMetadata[col][[1]]/sum(bw.norm@elementMetadata[col][[1]])
    # bw.norm@elementMetadata[col][[1]] <- (bw.norm@elementMetadata[col][[1]] - mean(bw.norm@elementMetadata[col][[1]]))/sd(bw.norm@elementMetadata[col][[1]])
    # bw.norm@elementMetadata[col][[1]] <- bw.norm@elementMetadata[col][[1]] +2
    # bw.norm@elementMetadata[col][[1]] <- (bw.norm@elementMetadata[col][[1]] - min(bw.norm@elementMetadata[col][[1]]))/(max(bw.norm@elementMetadata[col][[1]])-min(bw.norm@elementMetadata[col][[1]]))
      }
  return(bw.norm)
}

## Trying to plot
plot.expVSpred.coverage.track <- function(common.bw, met=NULL, plot=TRUE){
  gtrack <- GenomeAxisTrack()
  chrom <- as.character(unique(common.bw@seqnames))
  start <- common.bw@ranges@start[1]
  end <- common.bw@ranges@start[length(common.bw@ranges)]
  dtrack <- DataTrack(common.bw, 
                      chromosome=chrom,
                      name='norm. coverage',
                      groups=colnames(values(common.bw)[sapply(values(common.bw), is.numeric)]),
                      type='l'
                      )
  tracklist <- list(gtrack,dtrack)
  if (!is.null(met)) {
    metrack <- DataTrack(met,
                         chromosome=chrom,
                         name = 'methylation frac.',
                         type='gradient',
                         gradient=RColorBrewer::brewer.pal(5,'RdBu'),
                         ylim=c(0,100))
    tracklist <- list(gtrack,dtrack, metrack)
  }
  if (plot){
    plotTracks(tracklist, 
               from=start, to=end,
               # type=c('h', 'l', ), 
               main=paste0(chrom,':',start,':', end),
               # groups=colnames(values(common.bw)[sapply(values(common.bw), is.numeric)]), 
               legend=TRUE)
  }
  return(tracklist)
}

plot.subset <-function(bw,chr,start = 0,end = 1000000000000, met=NULL){
  subset.bw <- subsetByRegion(bw, chr,start,end)
  plot.expVSpred.coverage.track(subset.bw, met=met)
}

### Comparing peak heights

computeDistance <- function(bw.ranges){
  sim <- similarity(bw.ranges$score, bw.ranges$pred)
  met.df <- data.frame(spearman = sim$metrics$SPEARMAN_CORRELATION, 
             maxmax=sim$metrics$RATIO_INTERSECT, 
             id=unique(bw.ranges$id))
  return(met.df)
}

split.region <- function(gr, bps=100){
  sp <- split(gr, rep(1:round(length(gr)/bps),each=bps)[1:length(gr)])
  return(sp)
}

smooth.coverage <- function(gr, bandwith=100){
  gr$score <- ksmooth(1:length(gr$score), gr$score, kernel='normal', bandwidth = 100)$y
  return(gr)
}

adjust.prediction <- function(test.bw, plot=FALSE, lm.summary=FALSE){
  perc.rank <- function(x) trunc(rank(x))/length(x)
  local.norm.test.bw <- normalize.coverage(test.bw)
  df <- data.frame(values(local.norm.test.bw)) %>% arrange(-score) 
  df <- df %>% mutate(quant.pred = 1-perc.rank(pred)) %>% mutate(quant.score = 1-perc.rank(score))
  df$quant.pred[df$quant.pred==0] <- 1e-10
  mod.loc <- lm(score ~ log(quant.pred), data=df)
  s <- summary(mod.loc)
  if (lm.summary) {
    print(s)
  }
  ranked.values <- data.frame(values(local.norm.test.bw)) %>% 
    mutate(quant.pred = 1-perc.rank(pred))  %>% 
    mutate(quant.score = 1-perc.rank(score))
  ranked.values$quant.pred[ranked.values$quant.pred==0] <- NA
  values.model <- ranked.values %>% 
    mutate(adj.pred = predict(mod.loc, data.frame(quant.pred = ranked.values$quant.pred)))
  if (plot) {
    cols <- c('real.cov'='black', 'pred.cov'='blue', 'adj.pred.cov'='red')
    p <- ggplot(values.model, aes(quant.score, score)) +
      geom_point(aes(color='real.cov')) +
      geom_point(aes(quant.pred, pred, color='pred.cov')) +
      geom_line(aes(quant.pred, (adj.pred), color='adj.pred.cov')) +
      xlab('quantile') + 
      scale_color_manual(name='', values = cols) +
      annotate('text', x=0.9, y=max(values.model$score)- (5*max(values.model$score))/100, label=paste('R.sq =',round(s$adj.r.squared,2), '\n'), size=5) +
      # scale_y_log10() +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=10), title = element_text(size=22)) 
    print(p)
  }
  local.norm.test.bw$adj.pred <- values.model$adj.pred
  return(local.norm.test.bw)
}


adjust.prediction.exp <- function(test.bw, plot=F){
  perc.rank <- function(x) trunc(rank(x))/length(x)
  df <- as.data.frame(values(test.bw)) %>% 
    mutate(quant.score = perc.rank(x = score)) %>%
    mutate(quant.pred = perc.rank(x = pred)) 
  lm2 <- lm(log(score) ~ quant.score, 
            data=subset(df, df$score!=0), # To avoid error from -Inf
            na.action = 'na.exclude')
  if (plot) {
    p <- ggplot(df, aes(quant.score, score)) +
      geom_point() + 
      geom_point(color='red', aes(quant.pred, pred)) +
      geom_line(aes(quant.pred, exp(predict(lm2, data.frame(quant.score=quant.pred)))), color='red')
    
    print(p)
  }
  values(test.bw) <- df %>% 
    mutate(adj.pred =exp(predict(lm2, data.frame(quant.score=quant.pred)))) %>%
    select(score,pred,adj.pred, id)
  return(test.bw)
  }

get.range.methylation <- function(test.bw, baseres.met){
  ovs <- findOverlaps(test.bw, baseres.met)
  meth.prof <- baseres.met[subjectHits(ovs)]
  values(meth.prof) <- values(meth.prof)$frac
  # values(test.bw) <- data.frame(values(test.bw)) %>% mutate( met.frac=NA)
  # values(test.bw)$met.frac[queryHits(ovs)] <- meth.prof$frac
  return(meth.prof)
}

get.avg.methylation <- function(test.bw, baseres.met){
  met <- get.range.methylation(test.bw, baseres.met = baseres.met)$X
  return(mean(met))
}

plot.cov.wAnnotation <- function(test.bw, anno.gr){
  test.ovs <- findOverlaps(test.bw, anno.gr)
  p <- plot.expVSpred.coverage.track(normalize.coverage(test.bw), plot = FALSE)
  anno.track <- AnnotationTrack(anno.gr[subjectHits(test.ovs)],
                                name='CTCF')
  plotTracks(list(p[[1]], p[[2]], anno.track),
             main=unique(test.bw$id))
}

