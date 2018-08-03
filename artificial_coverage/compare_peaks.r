### Compare experimental and predicted coverage
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(nortest)
library(pheatmap)
library(rtracklayer)
# source("https://bioconductor.org/biocLite.R")
library(parallel)
library(GenomicRanges)

make.base.res.bw <- function(bw){
  base.res.bw.list <- tile(bw, width=1)
  base.res.bw <- base.res.bw.list@unlistData
  base.res.bw$score <- rep(bw$score,width(bw))
  return(base.res.bw)
}

make.predVSexp.range <- function(pred.bw, exp.bw, pred.name='pred', exp.name='score'){
  ranges <- subsetByOverlaps(exp.bw,pred.bw)
  hits <- findOverlaps(exp.bw, pred.bw)
  idx <- subjectHits(hits)
  values <- DataFrame(pred = pred.bw$score[idx])
  mcols(ranges) <- c(mcols(ranges), values)
  colnames(mcols(ranges)) <- c(exp.name, pred.name)
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

add.id.2 <- function(common.ranges, reg.length=3000){
  common.ranges$range.id <- cut(seq_along(common.ranges), breaks = seq(0, length(common.ranges), by = reg.length))
  return(common.ranges)
  }

test.add.id <- function(range, reg.length=3000){
  l <- max(range@ranges@start) - min(range@ranges@start)
  chroms <- length(range@seqnames@values)
  if(l==reg.length-1 & chroms==1){
    return('good')
  } else {
    return('bad')
    }
  }

split.tracks <- function(common.ranges, reg.length){
  range.list <- split(common.ranges, common.ranges$range.id)
  all <- sapply(range.list, test.add.id, reg.length=reg.length)
  if (any(all=='bad')){stop('Ranges are not all of the same length!')
  } else {
      return(range.list)
    }
  }


load.expVSpred.coverage <- function(pred.bw.file, exp.bw.file, save=FALSE){
  print('Loading predicted coverage profile...')
  pred.bw <- import(pred.bw.file, format = 'BigWig')
  ranges <- GRanges(seqnames = pred.bw@seqnames, 
                    ranges = IRanges(start=pred.bw@ranges@start, 
                                     end = pred.bw@ranges@start+1))
  print('Loading expected coverage profile...')
  exp.bw <- import(exp.bw.file, format = 'BigWig', which = ranges)
  exp.bw.base <- make.base.res.bw(exp.bw)
  # print('Merging coverage profiles...')
  # common.bw <- make.predVSexp.range(pred.bw = pred.bw, exp.bw = exp.bw.base)
  # common.bw <- add.id(common.bw)
  # if (save) {
  #   save(human_noBS.common.bw, file='~/AvOwork/human_noBS.covprofiles.RData')
  # }
  return(list(pred=pred.bw, exp=exp.bw))
  }

make.predVSexp.track <- function(predicted.bw, experimental.bw, reg.length=3000){
  pred.exp.tracks <- load.expVSpred.coverage(pred.bw.file = predicted.bw, exp.bw.file = experimental.bw)
  br <- smooth.coverage(make.base.res.bw(pred.exp.tracks[[2]]))
  track <- make.predVSexp.range(pred.exp.tracks[[1]], br) %>%
    normalize.coverage() %>%
    add.id.2(reg.length = reg.length)
  track.list <- lapply(split.tracks(track, reg.length = reg.length), trim.edges)
  return(track.list)  
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
    # bw.norm@elementMetadata[col][[1]] <- bw.norm@elementMetadata[col][[1]]/sum(bw.norm@elementMetadata[col][[1]])
    bw.norm@elementMetadata[col][[1]] <- (bw.norm@elementMetadata[col][[1]] - mean(bw.norm@elementMetadata[col][[1]], na.rm=TRUE))/sd(bw.norm@elementMetadata[col][[1]], na.rm=TRUE)
    # bw.norm@elementMetadata[col][[1]] <- bw.norm@elementMetadata[col][[1]] +2
    # bw.norm@elementMetadata[col][[1]] <- (bw.norm@elementMetadata[col][[1]] - min(bw.norm@elementMetadata[col][[1]]))/(max(bw.norm@elementMetadata[col][[1]])-min(bw.norm@elementMetadata[col][[1]]))
      }
  return(bw.norm)
}

## Trying to plot
plot.expVSpred.coverage.track <- function(common.bw, met=NULL, plot=TRUE){
  options(ucscChromosomeNames=FALSE)
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

# split.region <- function(gr, bps=100){
#   sp <- split(gr, rep(1:round(length(gr)/bps),each=bps)[1:length(gr)])
#   return(sp)
# }

smooth.coverage <- function(gr, bandwith=100){
  gr$score <- ksmooth(1:length(gr$score), gr$score, kernel='normal', bandwidth = 100)$y
  return(gr)
}

trim.edges <- function(gr, trim=100){
  trimmed.gr <- gr[(trim+1):(length(gr)-trim)]
  return(trimmed.gr)
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

nice.plotTrack <- function(test.bw, labels=c('Experimental ', 'Predicted')){
  chrom <- test.bw@seqnames@values
  long.df <- as.data.frame(values(test.bw)) %>% 
    select(-range.id) %>%
    mutate(genomic.coord=test.bw@ranges@start) %>%
    melt(id.vars='genomic.coord', variable.name='sample', value.name='norm.coverage')
  p <- ggplot(long.df, aes(genomic.coord, norm.coverage, group=sample, color=sample, linetype=sample)) +
    geom_line(size=2) +
    xlab("Genomic coordinates") +
    ylab('Norm. coverage (zscore)') +
    ggtitle(paste0(chrom, ':', as.character(min(long.df$genomic.coord)),":", as.character(max(long.df$genomic.coord)))) +
    theme_bw() +
    scale_color_manual(values=c('royalblue3', 'gold2'),labels=labels) +
    scale_linetype_manual(values=c(8,1), labels=labels) +
    theme(axis.title = element_text(size = 28), 
          axis.text = element_text(size=18), 
          plot.title = element_text(size=30, hjust=0.5), 
          legend.text=element_text(size=25),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.key.size = grid::unit(5, "lines"))
  return(p)
}

#### CORRELATION ANALYSIS ####
correlate.pred.exp <- function(pred, exp, method='spearman'){
  cor.coef <- cor(exp, pred, method=method)
  return(cor.coef)  
}

corr.real <- function(track, method='spearman'){
  return(correlate.pred.exp(exp=track$score, pred=track$pred)) 
}

randomize <- function(track.list.clean, method='spearman'){
  rand <- sample(track.list.clean)
  corr.coef.rand <- map_dbl(seq_along(track.list.clean), function(i) correlate.pred.exp(exp=track.list.clean[[i]]$score, pred=rand[[i]]$pred, method = method))
  return(corr.coef.rand)
}

compare.spear.real <- function(track.list, name.real='WGS', method='spearman'){
  spear <- map_dbl(track.list, corr.real,method=method)
  rand <- randomize(track.list, method=method)
  pl <- data.frame(sample=name.real, spearman.rho=spear) %>%
    bind_rows(., data.frame(sample='random', spearman.rho=rand)) %>%
    ggplot(., aes(sample, spearman.rho)) +
    geom_boxplot(outlier.alpha = 0.2, varwidth = T) +
    theme_classic() +
    ggsignif::geom_signif(comparisons = list(c(name.real, "random")),
                          test='wilcox.test',
                          map_signif_level=TRUE) +
    theme(axis.title=element_text(size=20), 
          axis.text=element_text(size=16), 
          strip.text=element_text(size=20), 
          legend.text = element_text(size=20),
          legend.title = element_text(size=22)) +
    xlab('')
  return(list(spear=spear, rand=rand, p=pl))
}
  

#### AUC YIELD ANALYSIS ####


coverage.yield.delta <- function(scaled.track, roi.track){
  int.roi <- findOverlaps(query = scaled.track, subject = roi.track)
  auc.best.roi <- auc(seq_along(scaled.track[queryHits(int.roi)]$best), scaled.track[queryHits(int.roi)]$best)
  auc.even.roi <- auc(seq_along(scaled.track[queryHits(int.roi)]$even), scaled.track[queryHits(int.roi)]$even)
  random.out <- scaled.track[-queryHits(int.roi)]
  auc.even.out <- auc(seq_along(sample(random.out$even, length(queryHits(int.roi)))), sample(random.out$even, length(queryHits(int.roi))))
  auc.best.out <- auc(seq_along(sample(random.out$best, length(queryHits(int.roi)))), sample(random.out$best, length(queryHits(int.roi))))
  yield.best <- auc.best.roi/auc.best.out
  yield.even <- auc.even.roi/auc.even.out
  return(yield.best-yield.even)
}

coverage.yield.single <- function(scaled.track, roi.track){
  int.roi <- findOverlaps(query = scaled.track, subject = roi.track)
  auc.best.roi <- auc(seq_along(scaled.track[queryHits(int.roi)]$score), scaled.track[queryHits(int.roi)]$score)
  random.out <- scaled.track[-queryHits(int.roi)]
  auc.best.out <- auc(seq_along(sample(random.out$score, length(queryHits(int.roi)))), sample(random.out$score, length(queryHits(int.roi))))
  yield.best <- auc.best.roi/auc.best.out
  return(yield.best)
}

norm.scale.n.yield <- function(raw.track, roi.track, scale.by=5){
  norm.track <- normalize.coverage(raw.track) %>%
    add.id.2(reg.length = 2000)
  scaled.norm.track <- norm.track
  score.cols <- colnames(values(scaled.norm.track)[sapply(values(scaled.norm.track), is.numeric)])
  for (col in score.cols) {
    scaled.norm.track@elementMetadata[col][[1]] <- scaled.norm.track@elementMetadata[col][[1]] + scale.by
  }
  yield <- coverage.yield.single(scaled.norm.track, roi.track)
  return(yield)
}

make.cum.dist <- function(vec, plot=T){
  cum.dist.best.out <- data.frame(even=vec) %>%
    mutate(cumdist=cume_dist(even)) 
  if (plot) {
    p <- cum.dist.best.out %>%
      ggplot(., aes(even, cumdist)) +
      geom_point(size=0.5)
  }
  return(list(cumdist=cum.dist.best.out, p=p))
}

random.from.even.2 <- function(vec, vec.cum.dist){
  r <- runif(1,0,1)
  s <- suppressWarnings(max(vec[r >= vec.cum.dist], na.rm=T))
  return(s)
}

make.random.profile <- function(cum.dist.real, threads = detectCores()){
  # random.profile <- c()
  # n <- 1
  random.profile <- mclapply(1:nrow(cum.dist.real$cumdist), function(x)
      random.from.even.2(cum.dist.real$cumdist$even, cum.dist.real$cumdist$cumdist),
      mc.cores = threads
  )
  # while (n <= nrow(cum.dist.real$cumdist)){
  #   r <- random.from.even.2(cum.dist.real$cumdist$even, cum.dist.real$cumdist$cumdist)
  #   random.profile <- c(random.profile, r )
  #   n <- n+1
  # }
  return(unlist(random.profile))
}

# bench.mark.random.profile <- function(cum.dist){
#   a <- sapply(1:7, function(n) system.time(p <- make.random.profile(cum.dist, threads=n))
#   }


delta.yield.permutation <- function(my.track, roi.track, permute='best', threads=detectCores()){
  rand.prof <- make.random.profile(make.cum.dist(my.track$even), threads = threads)
  rand.prof[is.infinite(rand.prof)] <- 0
  permuted.track <- my.track
  if (permute=='best') {
    permuted.track@elementMetadata$best <- rand.prof
  } else if (permute=='even'){
    permuted.track@elementMetadata$even <- rand.prof
  }
  yield.delta <- coverage.yield.delta(permuted.track, roi.track)
  return(yield.delta)  
}

random.delta.yield.dist <- function(my.track, roi.track, n.iterations=1000, threads=detectCores(), verbose=T){
  real.delta <- coverage.yield.delta(my.track, roi.track)
  it <- 1
  random.deltas <- c()
  while (it< n.iterations) {
    if (verbose) {  print(paste("Running iteration no.", it), quote = F)  }
    d <- delta.yield.permutation(my.track, roi.track, threads = threads)  
    random.deltas <- c(random.deltas, d)
    it <- it+1  
  }
  return(list(real=real.delta, random=random.deltas))
}
