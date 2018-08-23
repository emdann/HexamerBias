### Keq COMPUTATION ###
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
source('~/HexamerBias/rscripts/sanity_checks_deltaG.r')

# Putting together functions used in normalization_strategy.Rmd

loadPtMatrix <- function(file, compression='gzip'){
  if (compression=='gzip') {
    tab <- fread(paste('zcat',file), sep = ',', header=TRUE)
  }else{
    tab <- fread(file, sep = ',', header=TRUE) 
  }
  p <- as.matrix(tab[,-1])
  rownames(p)<- tab$V1
  return(p)
}

make.hex.usage.df <- function(ptTab, type = 'primer', scale =TRUE){
  if (type=='primer') { hex.usage <- colSums(ptTab)  }
  if (type=='template') { hex.usage <- rowSums(ptTab)  }
  sample.name <- gsub(deparse(substitute(ptTab)), pattern = 'pt.', replacement = '')
  if (scale) { 
    hex.df <- data.frame(names(hex.usage), hex.usage/sum(hex.usage), row.names = NULL) 
  }else{
    hex.df <- data.frame(names(hex.usage), hex.usage, row.names = NULL)
  }
  colnames(hex.df) <- c(type,sample.name)
  return(hex.df)
}


load.pt.data <- function(ptCounts.file, diag.pairs = F){
  # Loading primer-template matrix, returns long matrix 
  pt <- loadPtMatrix(ptCounts.file,compression = 'none')
  printf("primer-template matrix loaded!")
  template.usage <- make.hex.usage.df(pt, type='template', scale=F)
  print("template usage dataframe done!")
  primer.usage <- make.hex.usage.df(pt, type='primer', scale=F)
  print("primer usage dataframe done!")
  colnames(template.usage) <- c('template', 't.usage')
  colnames(primer.usage) <- c('primer', 'p.usage')
  pairs <- make_pair_df(pt)
  print("Long dataframe for pairs done!")
  if (diag.pairs) {
    matches <- pairs[substr(pairs$ptPair,1,6)==substr(pairs$ptPair,8,13),]
  } else {
    matches <- pairs
  }
  return(list(matches=matches, t.usage=template.usage, p.usage=primer.usage))
}

load.kmer.abundance <- function(kmer.abundance.file){
  return(read.csv(kmer.abundance.file, header = F, col.names = c('template', 'abundance')))
}

load.modelled.deltaG <- function(deltaG.file){
  return(read.csv(deltaG.file, header = F, col.names = c('template', 'dG')))
}

join.pt.data <- function(matches, template.usage, kmer.ab, tabDg){
  df <- matches %>%
    transmute(template=substr(ptPair,1,6), primer=substr(ptPair,8,100), pt=dG) %>%
    inner_join(., df, by='template') %>%
    inner_join(., template.usage, by='template')
  # rename(t.usage=pt)
  return(df)
}

select.diag.pairs <- function(pairs.df){
  pairs.df <- filter(pairs.df, template == primer)
  return(pairs.df)
}

pt.data <- load.pt.data(ptCounts.file = "~/mnt/edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-CeleTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv", diag.pairs = F)
