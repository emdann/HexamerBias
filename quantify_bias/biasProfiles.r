## Plot profiles
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)

get.profile.from.matrix <- function(file.gz){
  mat <- read.table(gzfile(file.gz), skip = 1)
  profMat <- t(mat[,-c(1:6)])
  profile <- apply(profMat,1,mean, na.rm=TRUE)
  return(profile/sum(profile))
}

VAN1667.profile <-profile/sum(profile)
VAN1667.subset <-get.profile.from.matrix("~/mnt/edann/hexamers/strand_specific/VAN1667_se.random42.srt.mat.gz")
profiles <- list(VAN1667.profile, VAN1667.subset)
noPBAT.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/kaester/ERR454965_1_val_1_bismark_bt2.deduplicated.srt.mat.gz")

make.df.of.profiles <- function(profiles){
  prof.df <- data.frame(profiles) %>%
    mutate(position = seq(1,1100)) %>%
    melt(id.vars = c('position'), variable.name='sample')
  return(prof.df)
}
# VAN1667.subset.profile.df <- data.frame(position = seq(1,1100),VAN1667.subset = VAN1667.subset)
profiles <- list(all=VAN1667.profile, subset=VAN1667.subset)
prof.df <- make.df.of.profiles(profiles)

ggplot(prof.df, aes(position,value, color=sample)) + 
  theme_classic() +
  geom_line() +
  scale_x_continuous(breaks = c(0,300,900,1100), 
                     labels = c('0' = '-3kb', '300' = 'TSS', '900' = 'TES', '1100' = '+3kb')) +
  xlab('Relative position') +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=10), title = element_text(size=22)) 
  
