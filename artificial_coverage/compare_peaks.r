###################################
### PREDICTED COVERAGE PROFILES ###
###################################
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

#' Genomic ranges to base resolution
#' 
#' Breaks GRanges object to single base resolution, giving same score to every base from a same region
#' 
#' @param bw GRanges object 
#' 
#' @return GRanges object of base resolution track
#' 
#' @export
make.base.res.bw <- function(bw){
  base.res.bw.list <- tile(bw, width=1)
  base.res.bw <- base.res.bw.list@unlistData
  base.res.bw$score <- rep(bw$score,width(bw))
  return(base.res.bw)
}

#' Merge coverage profiles
#' 
#' Generates a single GRanges object containing overlapping ranges for 2 GRanges objects with score (coverage)
#' 
#' @param pred.bw GRanges track (predicted coverage profile)
#' @param exp.bw GRanges track (experimental coverage profile)
#' @param pred.name character of name for first track score column (e.g. "predicted")
#' @param exp.name character of name for second track score column
#' 
#' @return GRanges track of common ranges with 2 score columns
#' 
#' @export
make.predVSexp.range <- function(pred.bw, exp.bw, pred.name='pred', exp.name='score'){
  ranges <- subsetByOverlaps(exp.bw,pred.bw)
  hits <- findOverlaps(exp.bw, pred.bw)
  idx <- subjectHits(hits)
  values <- DataFrame(pred = pred.bw$score[idx])
  mcols(ranges) <- c(mcols(ranges), values)
  colnames(mcols(ranges)) <- c(exp.name, pred.name)
  return(ranges)
} 

#' Add region ID
#' 
#' Deprecated. Use \code{add.id.2}
add.id <- function(ranges){
  breaks <- sapply(seq_along(ranges@ranges@start), function(ix)
    ranges@ranges@start[ix+1] - ranges@ranges@start[ix]!=1)
  breaks <- unlist(breaks)
  break.points <- unique(ranges[which(breaks)+1]@ranges@start)
  ids <- cut(ranges@ranges@start, breaks = c(0,break.points, max(break.points)+10000), include.lowest = TRUE)
  ranges$id <- paste(ranges@seqnames,ids)
  return(ranges)
  }


#' Add region ID
#' 
#' Adds column containing region id (chr:start:end) to every track line in base resolution GRanges object,
#' so that adjacent positions in the same region get the same ID. 
#' All ranges must have the same length.
#' 
#' @param ranges GRanges object
#' @param reg.length length of regions 
#' 
#' @return GRanges object with ID column
#' 
#' @export 
add.id.2 <- function(ranges, reg.length=3000){
  ranges$range.id <- cut(seq_along(ranges), breaks = seq(0, length(ranges), by = reg.length))
  return(ranges)
  }

#' Test IDs
#' 
#' Tests if region IDs correspond to start and end of region
#' 
#' @param range GRanges object with id column
#' @param reg.length length of regions
#' 
#' @return 'good' or 'bad'
#' 
#' @export
test.add.id <- function(range, reg.length=3000){
  l <- max(range@ranges@start) - min(range@ranges@start)
  chroms <- length(range@seqnames@values)
  if(l==reg.length-1 & chroms==1){
    return('good')
  } else {
    return('bad')
    }
  }

#' Split ranges by ID
#' 
#' Splits GRanges object in GRanges list of ranges with the same ID 
#' 
#' @param range GRanges object
#' @param reg.length length of regions
#' 
#' @return GRangesList object as long as the number of regions
#' 
#' @export
split.tracks <- function(range, reg.length){
  range.list <- split(range, range$range.id)
  all <- sapply(range.list, test.add.id, reg.length=reg.length)
  if (any(all=='bad')){stop('Ranges are not all of the same length!')
  } else {
      return(range.list)
    }
  }

#' Standard normalize coverage 
#' 
#' Compute z-score normalization (x-mean(x))/sd(x)) for coverage profiles 
#' 
#' @param bw GRanges object with score columns 
#' 
#' @return GRanges object with normalized coverage profiles 
#' 
#' @export
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

#' Load coverage profiles
#' 
#' Loads experimental and predicted coverage profiles from bigWig files 
#' 
#' @param pred.bw.file path to bigWig file of predicted coverage profile
#' @param exp.bw.file path to bigWig file of experimental coverage profile
#' 
#' @return list containing 2 GRanges objects
#' 
#' @export
load.expVSpred.coverage <- function(pred.bw.file, exp.bw.file){
  print('Loading predicted coverage profile...')
  pred.bw <- import(pred.bw.file, format = 'BigWig')
  ranges <- GRanges(seqnames = pred.bw@seqnames, 
                    ranges = IRanges(start=pred.bw@ranges@start, 
                                     end = pred.bw@ranges@start+1))
  print('Loading expected coverage profile...')
  exp.bw <- import(exp.bw.file, format = 'BigWig', which = ranges)
  exp.bw.base <- make.base.res.bw(exp.bw)
  return(list(pred=pred.bw, exp=exp.bw))
  }

#' Kernel smoothing of coverage profiles 
#' 
#' @param gr GRanges object 
#' @param bandwidth the bandwidth. The kernels are scaled so that their quartiles (viewed as probability densities) are at +/- 0.25*bandwidth.
#' 
#' @return GRanges object of smoothened coverage profile
#' 
#' @export
smooth.coverage <- function(gr, bandwith=100){
  gr$score <- ksmooth(1:length(gr$score), gr$score, kernel='normal', bandwidth = 100)$y
  return(gr)
}

trim.edges <- function(gr, trim=100){
  trimmed.gr <- gr[(trim+1):(length(gr)-trim)]
  return(trimmed.gr)
}


#' Make track of predicted VS experimental coverage 
#' 
#' Merges, smoothens and normalizes coverage tracks for experimental and predicted coverage, adding IDs for each sampled region
#' (region length has to be specified)
#' 
#' @param predicted.bw GRanges object of predicted coverage track
#' @param experimental.bw GRanges object of experimental coverage track
#' @param reg.length length of sampled regions 
#' 
#' @return GRangesList object where each element is a sampled region
#' 
#' @export
make.predVSexp.track <- function(predicted.bw, experimental.bw, reg.length=3000){
  pred.exp.tracks <- load.expVSpred.coverage(pred.bw.file = predicted.bw, exp.bw.file = experimental.bw)
  br <- smooth.coverage(make.base.res.bw(pred.exp.tracks[[2]]))
  track <- make.predVSexp.range(pred.exp.tracks[[1]], br) %>%
    normalize.coverage() %>%
    add.id.2(reg.length = reg.length)
  track.list <- lapply(split.tracks(track, reg.length = reg.length), trim.edges)
  return(track.list)  
}

### Subset
#' Subset by region
#' 
#' @param ranges GRanges object 
#' @param chrom selected chromosome 
#' @param start selected starting position
#' @param end selected ending position
#' 
#' @return GRanges object of track in selected coordinates
#' 
#' @export
subsetByRegion <- function(ranges, chrom, start, end){
  subset.bw <- ranges[ranges@seqnames==chrom &
                  ranges@ranges@start >= start &
                  ranges@ranges@start < end]
  return(subset.bw)  
}

## Trying to plot

#' Plot predicted and experimental coverage profiles 
#' 
#' Deprecated, use niceplotTrack
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

#' Extract methylation of genomic region
#' 
#' @param reg.track GRanges track of one sampled region (single element of GRanges List slit by ID)
#' @param baseres.met GRanges object of base resolution methylation fraction
#' 
#' @return GRanges object of methylation values over the input region
#' 
#' @export
get.range.methylation <- function(reg.track, baseres.met){
  ovs <- findOverlaps(reg.track, baseres.met)
  meth.prof <- baseres.met[subjectHits(ovs)]
  values(meth.prof) <- values(meth.prof)$frac
  return(meth.prof)
}

#' Compute average methylation of genomic region
#' 
#' @param reg.track GRanges track of one sampled region (single element of GRanges List slit by ID)
#' @param baseres.met GRanges object of base resolution methylation fraction
#' 
#' @return average meth fraction in input region 
#' 
#' @export
get.avg.methylation <- function(reg.track, baseres.met){
  met <- get.range.methylation(reg.track, baseres.met = baseres.met)$X
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

#' Plot comparison of coverage profiles
#' 
#' @description Make ggplot viz of predicted and experimental coverage tracks for a sampled region of the genome 
#' 
#' @param region.track GRanges track of one sampled region (single element of GRanges List slit by ID)
#' @param labels vector of names for the 2 compared coverage profiles 
#' 
#' @return ggplot object 
#' 
#' @export
nice.plotTrack <- function(region.track, labels=c('Experimental ', 'Predicted')){
  chrom <- region.track@seqnames@values
  long.df <- as.data.frame(values(region.track)) %>% 
    select(-range.id) %>%
    mutate(genomic.coord=region.track@ranges@start) %>%
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
#' Correlation of predicted and experimental coverage 
#' 
#' @param pred vector of predicted coverage profile
#' @param exp vector of experimental coverage profile
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson", "kendall", or "spearman" (default)
#' 
#' @return correlation coefficient 
#' 
#' @export
correlate.pred.exp <- function(pred, exp, method='spearman'){
  cor.coef <- cor(exp, pred, method=method)
  return(cor.coef)  
}

#' Compute coverage correlation of all sampled regions
#' 
#' @param track.list GRangesList object of sampled regions with experimental and predicted coverage 
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson", "kendall", or "spearman" (default)
#' 
#' @return vector of correlation coefficients for all sampled regions
corr.real <- function(track.list, method='spearman'){
  return(correlate.pred.exp(exp=track.list$score, pred=track.list$pred)) 
}

#' Compute coverage correlation of random permutations of predicted coverage
#' 
#' @description Given a GRangesList of sampled regions, where for each sample the experimental and predicted coverage profile 
#' @description is given, the function shuffles the predicted coverage profiles, assigning them to random regions and computes the correlation coefficient
#' 
#' @param track.list GRangesList object of sampled regions with experimental and predicted coverage 
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson", "kendall", or "spearman" (default)
#' 
#' @return vector of correlation coefficients of permutations of sampled regions
#'
#' @export
randomize <- function(track.list, method='spearman'){
  rand <- sample(track.list)
  corr.coef.rand <- map_dbl(seq_along(track.list), function(i) correlate.pred.exp(exp=track.list[[i]]$score, pred=rand[[i]]$pred, method = method))
  return(corr.coef.rand)
  }

#' Quantify correlation of predicted and experimental coverage
#' 
#' For a given GRangesList, computes correlation between pred and exp profiles in every regions, in randomnly perputed regions 
#' and plots boxplot of comparisons, with p-value for Wilcoxon's test.
#' 
#' @param track.list GRangesList object of sampled regions with experimental and predicted coverage 
#' @param name.real character, name of seq library
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson", "kendall", or "spearman" (default)
#' 
#' @return list of vector for real correlation coefficients, for coefficients of random permutation, ggplot object of boxplot
#' 
#' @export
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

#' Compute difference in yield between batches
#' 
#' Computes difference in coverage yield for region of interest between predicted coverage profiles
#' of two different primer batches
#' 
#'  @param scaled.track GRanges object with values for predicted coverage of best batch (best) and random batch (even) (scaled, so no values < 0)
#'  @param roi.track GRanges object of regions of interest
#'  
#'  @return difference in yield
#'  
#'  @export
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

#' Compute coverage yield
#' 
#' Computes ratio between area under the curve in region of interest and outside of region of interest 
#' (sample of the same bps of the region of interest)
#' 
#' @param scaled.track GRanges object with values for predicted coverage of best batch (best) and random batch (even) (scaled, so no values < 0)
#' @param roi.track GRanges object of regions of interest
#' 
#' @return score for coverage yield
#' 
#' @export
coverage.yield.single <- function(scaled.track, roi.track){
  int.roi <- findOverlaps(query = scaled.track, subject = roi.track)
  auc.best.roi <- auc(seq_along(scaled.track[queryHits(int.roi)]$score), scaled.track[queryHits(int.roi)]$score)
  random.out <- scaled.track[-queryHits(int.roi)]
  auc.best.out <- auc(seq_along(sample(random.out$score, length(queryHits(int.roi)))), sample(random.out$score, length(queryHits(int.roi))))
  yield.best <- auc.best.roi/auc.best.out
  return(yield.best)
}

#' Normalize and compute yield 
#' 
#' Computes ratio between area under the curve in region of interest and outside of region of interest 
#' (sample of the same bps of the region of interest) after normalization and scaling to have no coverage values < 0
#' 
#' @param raw.track GRanges object with values for predicted coverage of best batch (best) and random batch (even)
#' @param roi.track GRanges object of regions of interest
#' @param scale.by scaling fator to add to each coverage value to have all > 0
#' 
#' @return score for coverage yield
#' 
#' @export
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

#' Build cumulative distribution of coverage profile 
#' 
#' To make random permutation of coverage profile
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

#' Make random permutation of coverage profile for random batch
random.from.even.2 <- function(vec, vec.cum.dist){
  r <- runif(1,0,1)
  s <- suppressWarnings(max(vec[r >= vec.cum.dist], na.rm=T))
  return(s)
}

#' Make random permutation of coverage profile
#' 
#' To calculate p-value of difference in yield
#' 
#' @param cum.dist.real cumulative distribution of coverage scores 
#' @param threads amound of threads to use
#' 
#' @return vector of permuted coverage profile
#' 
#' @export
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

#' Compute difference in yield of permuted coverage profile
#' 
#' @export
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

#' Random yield 
#' 
#' Computes difference in coverage yield for region of interest between predicted coverage profiles
#' of two different primer batches. To calculate p-value of yield score, it also computes the distribution of 
#' difference in yield between randomnly shuffled coverage profiles.
#' 
#' @param my.track GRanges object with values for predicted coverage of best batch (best) and random batch (even) (SCALED so no cov < 0)
#' @param roi.track GRanges object of regions of interest
#' @param n.iterations number of iterations 
#' @param threads number of threads to use
#' @param verbose logical indicating wheather to print progress
#' 
#' @return list of delta yield and vector of delta yield for random permutations 
#' 
#' @export
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
