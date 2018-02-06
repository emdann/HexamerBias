### SPREAD ANALYSIS ###
library(dplyr)
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("seqLogo")
library(ggplot2)
library(seqLogo)
require(ggseqlogo)

hex.df<-read.csv("~/hex_analysis.csv")

p <- hex.df %>% mutate(abundance_lowConv=abundance_lowEff$abundance[match(hex.df.hex, abundance_lowEff$hex)]) %>% mutate(log_amp=log(usage/abundance_lowConv))
cl <- hclust(dist(p[,c("deltaG", "log_amp")]), method = "average")
p %>% mutate(cluster=as.character(cutree(cl,k=10))) %>% ggplot(.,aes(deltaG,log_amp, color=cluster)) + geom_point()
p %>% mutate(model=predict(lm(log_amp ~ deltaG, data=.)))  %>% ggplot(.,aes(deltaG,log_amp)) + geom_point() + geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], col="red")
mod <- p %>% lm(log_amp ~ deltaG, data=.)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
p %>% mutate(class=ifelse(deltaG < (-int +log_amp)/slope,"under","over")) %>% ggplot(.,aes(deltaG,log_amp, colour=class)) + geom_point() 
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
  div_p <- p   %>% mutate(class=ifelse(deltaG < (-int+thresh + log_amp)/slope & deltaG > (-int+thresh+ 1 + log_amp)/slope,"over","under"))
  over_hexs <- filter(div_p, class=="over")$hex.df.hex
  return(multiplot(ggplot(div_p,aes(deltaG,log_amp, colour=class)) + geom_point() + theme(legend.position='none'), makeSeqLogo(over_hexs), cols = 2))
}

pdf("~/Van Oudenaarden code/output/spread_logos_loxEff.pdf", onefile=T, height = 3, width = 7)
for(t in seq(-8,5)){  seqLogo_thresh(t) }
dev.off()

GCcont <- function(string){
  letters <- unlist(strsplit(as.character(string),''))
  count <- table(factor(letters, levels=c("A","T","C","G")))
  gc_cont <- sum(count[c("G","C")])/length(letters)
  return(gc_cont)
}

GCratio<- function(string){
  letters <- unlist(strsplit(as.character(string),''))
  count <- table(factor(letters, levels=c("A","T","C","G")))
  gc_ratio <- (count["C"]/length(letters))/GCcont(string)
  return(gc_ratio)
}

kmers <- read.delim("~/mnt/edann/count_kmers_wBottom_mm10.txt", col.names = c("hex", "count"))
cbind(hex.df, unconv_abundance = kmers[match(hex.df$hex.df.hex,kmers$hex),"count"]) %>% mutate(log_amp=log(usage/abundance_amp)) %>% mutate(ratio_conv_nonconv=abundance_amp/unconv_abundance) %>% mutate(ratio=ifelse(ratio_conv_nonconv>50, NA, ratio_conv_nonconv)) %>% ggplot(., aes(x=deltaG, y=log_amp, color=ratio, label=X)) + geom_point() + scale_colour_gradient2(midpoint = 1) + scale_shape_manual(name = "Treatment & State", labels=">50", )

cbind(hex.df, unconv_abundance = kmers[match(hex.df$hex.df.hex,kmers$hex),"count"]) %>% mutate(log_amp=log(usage/abundance_amp)) %>% mutate(ratio_conv_nonconv=abundance_amp/unconv_abundance) %>% mutate(ratio=ifelse(ratio_conv_nonconv>50, NA, ratio_conv_nonconv)) %>% ggplot(., aes(x=deltaG, y=log_amp, fill=ratio, color="", label=X)) + geom_point() + scale_colour_gradient2(midpoint = 1) +              
  scale_colour_manual(values=NA) +
  guides(colour=guide_legend("No data", override.aes=list(colour="black")))

abundance_lowEff <- read.delim("~/mnt/edann/hexamers/rand_hex_abundance_lowEfficiency.txt", col.names = c("hex", "abundance"), header = F)
abundance_lowEff <- abundance_lowEff[match(hex.df$hex.df.hex, abundance_lowEff$hex),]

pdf("~/Van Oudenaarden code/output/occupancy_deltaG_lowEffConv.pdf")
hex.df %>% mutate(abundance_lowConv=abundance_lowEff$abundance[match(hex.df.hex, abundance_lowEff$hex)]) %>% 
  mutate(log_amp=log(usage/abundance_lowConv)) %>% 
  ggplot(.,aes(deltaG,log_amp, label=hex.df.hex)) + 
  geom_point() + 
  # geom_text(cex=2, nudge_y = 0.4) + 
  ggtitle("Low efficiency conversion") + ylab("log(usage/abundance)")
dev.off()

mod <- hex.df %>% mutate(abundance_lowConv=abundance_lowEff$abundance[match(hex.df.hex, abundance_lowEff$hex)]) %>% 
  mutate(log_amp=log(usage/abundance_lowConv)) %>% lm(log_amp ~ deltaG, data=.)

dens_p <- hex.df %>% mutate(abundance_lowConv=abundance_lowEff$abundance[match(hex.df.hex, abundance_lowEff$hex)]) %>% 
  mutate(log_amp=log(usage/abundance_lowConv)) %>% 
  ggplot(.,aes(deltaG,log_amp, label=hex.df.hex)) + 
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  geom_point(alpha = 0.1, shape = 20) +
  # geom_text(cex=2, nudge_y = 0.4) + 
  ggtitle("Low efficiency conversion") + ylab("log(usage/abundance)") + labs(fill="Density")

pdf("~/Van Oudenaarden code/output/occupancy_deltaG_lowEffConv_density.pdf")
dens_p
dev.off()

