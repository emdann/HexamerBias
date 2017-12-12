### SPREAD ANALYSIS ###
library(dplyr)
library(data.table)
# source("https://bioconductor.org/biocLite.R")
# biocLite("seqLogo")
library(ggplot2)
library(seqLogo)
require(ggseqlogo)

hex.df<-read.csv("~/hex_analysis.csv")

p <- hex.df %>% mutate(log_amp=log(usage/abundance_amp)) %>% filter(log_amp!=Inf)
cl <- hclust(dist(p[,c("deltaG", "log_amp")]), method = "average")
p %>% mutate(cluster=as.character(cutree(cl,k=10))) %>% ggplot(.,aes(deltaG,log_amp, color=cluster)) + geom_point()
p %>% mutate(model=predict(lm(log_amp ~ deltaG, data=.)))  %>% ggplot(.,aes(deltaG,log_amp)) + geom_point() + geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], col="red")
mod <- p %>% lm(log_amp ~ deltaG, data=.)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
p %>% mutate(model=predict(lm(log_amp ~ deltaG, data=.)))  %>% mutate(class=ifelse(deltaG < (-int-1 + log_amp)/slope,"under","over")) %>% ggplot(.,aes(deltaG,log_amp, colour=class)) + geom_point() 
div_p <- p %>% mutate(model=predict(lm(log_amp ~ deltaG, data=.)))  %>% mutate(class=ifelse(deltaG < (-int-1.5 + log_amp)/slope,"under","over"))

## Make sequence logos
over_hexs <- filter(div_p, class=="over")$hex.df.hex
under_hexs <- filter(div_p, class=="under")$hex.df.hex

computeSeqLogo <- function(hexs){
  counts <- lapply(c(1,2,3,4,5,6), function(i) sapply(hexs, function(x) substr(x,i,i)))
  pfm <- data.frame(do.call(cbind, lapply(lapply(counts,factor, levels=c("A","T","C","G")),table)))
  ppm <- pfm/colSums(pfm)
  pwm <- makePWM(ppm)
  return(seqLogo(pwm,xaxis = F))
}

makeSeqLogo <- function(hexs){
  p1 <- ggseqlogo(as.character(hexs))
  p2 <- p1 + ylim(c(0,2))
  return(p2)
  }


seqLogo_thresh <-function(thresh){
  div_p <- p %>% mutate(model=predict(lm(log_amp ~ deltaG, data=.)))  %>% mutate(class=ifelse(deltaG < (-int+thresh + log_amp)/slope & deltaG > (-int+thresh+ 1 + log_amp)/slope,"over","under"))
  over_hexs <- filter(div_p, class=="over")$hex.df.hex
  return(multiplot(ggplot(div_p,aes(deltaG,log_amp, colour=class)) + geom_point() + theme(legend.position='none'), makeSeqLogo(over_hexs), cols = 2))
}

pdf("~/Van Oudenaarden code/output/spread_logos.pdf", onefile=T, height = 3, width = 7)
for(t in seq(-5,4)){  seqLogo_thresh(t) }
dev.off()

GCcont <- function(string){
  letters <- unlist(strsplit(as.character(string),''))
  count <- table(factor(letters, levels=c("A","T","C","G")))
  gc_cont <- sum(count[c("G","C")])/length(letters)
  return(gc_cont)
  }

kmers <- read.delim("~/mnt/edann/count_kmers_wBottom_mm10.txt", col.names = c("hex", "count"))
cbind(hex.df, unconv_abundance = kmers[match(hex.df$hex.df.hex,kmers$hex),"count"]) %>% mutate(log_amp=log(usage/abundance_amp)) %>% mutate(ratio_conv_nonconv=abundance_amp/unconv_abundance) %>% mutate(ratio=ifelse(ratio_conv_nonconv>50, NA, ratio_conv_nonconv)) %>% ggplot(., aes(x=deltaG, y=log_amp, color=ratio, label=X)) + geom_point() + scale_colour_gradient2(midpoint = 1) + scale_shape_manual(name = "Treatment & State", labels=">50", )

cbind(hex.df, unconv_abundance = kmers[match(hex.df$hex.df.hex,kmers$hex),"count"]) %>% mutate(log_amp=log(usage/abundance_amp)) %>% mutate(ratio_conv_nonconv=abundance_amp/unconv_abundance) %>% mutate(ratio=ifelse(ratio_conv_nonconv>50, NA, ratio_conv_nonconv)) %>% ggplot(., aes(x=deltaG, y=log_amp, fill=ratio, color="", label=X)) + geom_point() + scale_colour_gradient2(midpoint = 1) +              
  scale_colour_manual(values=NA) +
  guides(colour=guide_legend("No data", override.aes=list(colour="black")))
