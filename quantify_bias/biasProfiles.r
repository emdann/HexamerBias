## Plot profiles
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)

get.profile.from.matrix <- function(file.gz){
  mat <- read.table(file.gz, skip = 1)
  profMat <- t(mat[,-c(1:6)])
  profile <- apply(profMat,1,mean, na.rm=TRUE)
  return(profile/sum(profile))
}

VAN1667.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/VAN1667prediction/artificial_coverage/matrixVAN1667.gz")
# VAN1667.subset <-get.profile.from.matrix("~/mnt/edann/hexamers/strand_specific/VAN1667_se.random42.srt.mat.gz")
noPBAT.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/kaester/ERR454965_1_val_1_bismark_bt2.deduplicated.srt.mat.gz")
# noPBAT.small.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/kaester/ERR454965_1_val_1_bismark_bt2.deduplicated.srt.small.mat.gz")

load.profile <- function(profile.txt){
  profile <- scan(profile.txt)
  return(profile/sum(profile))
  }

make.df.of.profiles <- function(profiles){
  prof.df <- data.frame(profiles) %>%
    mutate(position = seq(1,1100)) %>%
    melt(id.vars = c('position'), variable.name='sample')
  return(prof.df)
}
# VAN1667.subset.profile.df <- data.frame(position = seq(1,1100),VAN1667.subset = VAN1667.subset)
profiles <- list(PBAT=VAN1667.profile, no.PBAT=noPBAT.profile)
prof.df <- make.df.of.profiles(profiles)

plot.profile.df <- function(df){
  ggplot(df, aes(position,value, color=sample)) + 
    theme_classic() +
    geom_line() +
    scale_x_continuous(breaks = c(0,300,800,1100), 
                       labels = c('0' = '-3kb', '300' = 'TSS', '800' = 'TES', '1100' = '+3kb')) +
    xlab('Relative position') + ylab('normalized coverage') +
    theme(axis.title = element_text(size = 20), 
          axis.text = element_text(size=10), 
          title = element_text(size=22),
          legend.text = element_text(size=16)) 
}

### BS vs noBS
pbat.profile <- load.profile('~/mnt/edann/noPreAmp_crypts/PBAT_mm10.profile.txt')
noBS.profile <-load.profile('~/mnt/edann/noPreAmp_crypts/noBs_profile.txt')

df <- make.df.of.profiles(list(BS=pbat.profile, noBS = noBS.profile))

## Purified vs non-purified
purified.profile <- load.profile("~/mnt/edann/hexamers/OUD2086prediction/OUD2086_distal.profile.txt")
df <- make.df.of.profiles(list(purified=purified.profile))

## artificial coverage
artCov.prof <- get.profile.from.matrix('~/mnt/edann/hexamers/strand_specific/artificial_coverage/mm10.random.42.artCov.mat.gz')
load("~/AvOwork/PBAT_control_coverageprof_refGen.RData")
df <- make.df.of.profiles(list(artificial_coverage = artCov.prof))
p <- rbind(df,prof.df) %>% filter(sample!='no.PBAT')
