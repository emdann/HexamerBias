### Plot optimization performance
# install.packages('jsonlite')
# library(jsonlite)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)

reshape.prob.mat <- function(prob.mat){
  # nuc <- c("A", "T", "C", "G")
  # pos<-seq(1,6)
  # colnames(prob.mat) <- as.vector(t(sapply(seq(1,6), function(pos) paste0(nuc,'.',pos))))
  long.prob.mat <- prob.mat %>% mutate(iter=rownames(prob.mat)) %>%
    melt(variable.name = 'nuc.pos', value.name = 'prob') %>% 
    mutate(nuc=substr(nuc.pos,1,1), pos=substr(nuc.pos,3,3)) %>%
    select(iter,prob,nuc,pos)
  # long.prob.mat$nuc <- factor(long.prob.mat$nuc, levels=c('A', 'T', 'C', 'G'))
  return(long.prob.mat)  
}

plot.nuc.matrix <- function(prob.row){ # input is long df of one iteration
  cols <- brewer.pal(4,'Dark2')
  mat <- prob.row %>%
    mutate(nuc=ifelse(prob==0,NA,nuc)) %>% filter(!is.na(nuc))  # To avoid having label for zeroes
  mat$nuc <- factor(mat$nuc, levels=c('A', 'T', 'C', 'G'))
  p <- ggplot(mat, aes(pos,prob, fill=nuc)) + 
    geom_bar(stat = 'identity', alpha=0.4) +
    scale_fill_manual(name='Base', values = cols) +
    # geom_logo(aes(label=nuc, y=prob))
    geom_text(aes(label=nuc, color=nuc, size=prob), position = position_stack(vjust = 0.5), show.legend = FALSE, alpha=1) +
    scale_radius(range=c(1,22)) +
    scale_color_manual(values = cols) +
    xlab('Position') + ylab('Fraction') +
    ggtitle(paste('Iteration no.', prob.row$iter[1])) +
    theme_classic() +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size=10), title = element_text(size=22)) 
  return(p)  
  }

plot.iteration <- function(long.mat, it){
  iter.mat <- filter(long.mat, iter==as.character(it))  
  return(plot.nuc.matrix(iter.mat))
}

make.matrix.gif.frames <- function(prob.mat, path, name){
  long.mat <- reshape.prob.mat(prob.mat)
  l <- lapply(seq(1,10), function(i) plot.iteration(long.mat, i))
  for(el in l){
    png(paste0(path,name,"_it", el$data$iter[1], "_matrix.png"))
    plot(el)
    dev.off()
  }
}

make.rho.df <- function(opt.output){
  rho.df <- data.frame(iteration = seq_along(opt.output$score),rho = 1-opt.output$score)
  return(rho.df)
}


## Compare best & target



