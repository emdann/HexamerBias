#### OPTIMIZATION FUNCTIONS ####
library(parallel)
library(rtracklayer)
library(purrr)
library(Gviz)
library(zoo)
library(flux)
source("~/HexamerBias/deltaGprediction/binding_model_functions.r")

compute.enrichment.score <- function(dens.df, fc.df){
  suppressWarnings(
    enrichment.score <- dens.df %>%
    inner_join(., fc.df, by="template") %>%
    mutate(dens.fc=binding.dens*fc) %>%
    summarise(score=sum(dens.fc)) %>%
    .$score
   )
  return(enrichment.score)
}

density.combo <- function(prob.vec, keqs.df=d3r.keqs, eps=epsilon.d3r){
  b.probs <- batch.prob.uniform(hexs=all.hexamers(), nuc.probs = prob.vec)
  pred.cov.b <- predict.coverage(keqs.df, eps, prob = b.probs)
  dens.df <- pred.cov.b %>%
    mutate(binding.dens = pred.cov/abundance)  %>% 
    dplyr::select(template, binding.dens)
  colnames(dens.df)[2] <- do.call(paste,c('dens',prob.vec, sep='_'))
  return(dens.df)
}

joining.fun <- function(...){
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  col1 = colnames(df1)[1]
  col2 = colnames(df2)[1]
  xxx = left_join(..., by = setNames(col2,col1))
  return(xxx)
}


plot.ranked.fc <- function(cgi.fc){
  p <- cgi.fc %>%
    arrange(-fc) %>%
    mutate(rank=min_rank(fc)) %>%
    mutate(label=ifelse(rank<=5 | rank>=max(rank)-5, as.character(template), '')) %>%
    ggplot(., aes(rank, fc, label=label)) + 
    geom_point(size=0.5, alpha=0.2) +
    geom_text_repel()
  return(p)
}

plot.ranked.score <- function(cgi.fc){
  p <- cgi.fc %>%
    arrange(-score) %>%
    mutate(rank=min_rank(score)) %>%
    mutate(even=ifelse(template=='dens_0.25_0.25_0.25_0.25', score, NA)) %>%
    mutate(even.label=ifelse(template=='dens_0.25_0.25_0.25_0.25', 'Random', NA)) %>%
    mutate(label=ifelse(rank>=max(rank)-5, as.character(template), '')) %>%
    ggplot(., aes(rank, score, label=label)) + 
    geom_point(size=0.5, alpha=0.2) +
    geom_text_repel() +
    geom_point(aes(y=even), color='red', size=2) +
    geom_text_repel(aes(label=even.label), color='red', nudge_y = 2)
  return(p)
}

fc.scores <- function(dens.table, fc.df){
  dens.table.l <- map(column_to_rownames(dens.table, var='template'), function(x) data.frame(template=dens.table$template, binding.dens=x))
  scores <- map_dbl(dens.table.l, compute.enrichment.score, fc.df=fc.df)
  return(scores)
}
