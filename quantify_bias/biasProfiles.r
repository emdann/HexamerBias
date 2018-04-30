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

VAN1667.profile <- load.profile("~/mnt/edann/hexamers/strand_specific/VAN1667.profile.txt")
# # VAN1667.subset <-get.profile.from.matrix("~/mnt/edann/hexamers/strand_specific/VAN1667_se.random42.srt.mat.gz")
# noPBAT.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/kaester/ERR454965_1_val_1_bismark_bt2.deduplicated.srt.mat.gz")
# # noPBAT.small.profile <- get.profile.from.matrix("~/mnt/edann/hexamers/kaester/ERR454965_1_val_1_bismark_bt2.deduplicated.srt.small.mat.gz")

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

plot.profile.df <- function(df){
  p <-ggplot(df, aes(position,value, color=sample)) + 
    theme_classic() +
    geom_line(size=2) +
    scale_x_continuous(breaks = c(0,300,800,1100), 
                       labels = c('0' = '-3kb', '300' = 'TSS', '800' = 'TES', '1100' = '+3kb')) +
    xlab('Relative position') + ylab('normalized coverage') +
    theme(axis.title = element_text(size = 30), 
          axis.text = element_text(size=25), 
          title = element_text(size=22),
          legend.text = element_text(size=25), 
          legend.key.size = unit(1,"cm"))
  return(p)
}

### BS vs noBS
pbat.profile <- load.profile('~/mnt/edann/noPreAmp_crypts/PBAT_mm10.profile.txt')
noBS.profile <-load.profile('~/mnt/edann/noPreAmp_crypts/noBS_noChrUn.profile.txt')

df <- make.df.of.profiles(list(BS=pbat.profile, noBS = noBS.profile))
plot.profile.df(df)

## Purified vs non-purified
VAN1667.profile <- load.profile("~/mnt/edann/hexamers/strand_specific/VAN1667.profile.txt")
purified.profile <- load.profile("~/mnt/edann/hexamers/OUD2086prediction/10_R1.profile.txt")
df <- make.df.of.profiles(list(purified=purified.profile, non.purified=VAN1667.profile))
plot.profile.df(df)

## artificial coverage
artCov.prof <- load.profile('~/mnt/edann/hexamers/strand_specific/artificial_coverage/mm10.random.42.1000smp.artCov.profile.txt')
VAN1667.subsmp.prof <- load.profile('~/mnt/edann/hexamers/strand_specific/VAN1667_se.random42.srt.profile.txt')

df <- make.df.of.profiles(list(artificial_coverage = artCov.prof, VAN1667 = VAN1667.subsmp.prof))
plot.profile.df(df)
ggsave("~/AvOwork/output/artificial_coverage/bias_artCovVSVAN1667subsmp.pdf")
## Hand-mixed profiles
CP <- load.profile("~/mnt/edann/crypts_bs/VAN2408/CP.profile.txt")
MP <- load.profile("~/mnt/edann/crypts_bs/VAN2408/MP.profile.txt")
plot.profile.df(make.df.of.profiles(list(hand.mixed.CP = CP, hand.mixed.MP = MP, machine.mixed = VAN1667.profile)))
ggsave("~/AvOwork/output/coverage_bias/hadMixVSmachineMix_covprofile.pdf")

## Priming VS ligation
lig <- load.profile("~/mnt/edann/SRR1769256_chr1.profile.txt")
prim <- load.profile("~/mnt/edann/VAN1667.chr1.profile.txt")
p <- plot.profile.df(make.df.of.profiles(list(priming=prim, ligation=lig)))
