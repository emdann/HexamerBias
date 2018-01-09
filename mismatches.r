library(dplyr)
library(ggplot2)
al <- read.delim("~/mnt/edann/alignment_L1.chr10", header=F)

pdf("~/Van Oudenaarden code/output/align_hist_test.pdf")
hist(al$V2, xlab = "Alignment score", main="Alignment sequenced hexamer and primed region")
dev.off()

hex.df<-read.csv("~/hex_analysis.csv")
hex.df$hex.df.hex
mean_score <- function(hex){return(mean(al$V2[which(al$V1==hex)]))}
a=sapply(as.character(hex.df$hex.df.hex), mean_score)

pdf("~/Van Oudenaarden code/output/align_deltaG.pdf")
p <- hex.df %>% mutate(mean_aligment_score=a) %>% 
ggplot(., aes(x=deltaG, y=mean_aligment_score)) + geom_point()
p
dev.off()

hist <- read.csv("~/mnt/edann/bigHist_chr10_L1.csv", row.names = 1, header=T)
hist$N<-NULL
hist<- hist[-5,]
names(hist)[1]<-"-"
pheatmap(as.matrix(hist),cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 12, fontsize_col = 12,ylab="Base in sequenced hex", xlab="Base in primed region", 
         filename ="~/Van Oudenaarden code/output/mismatchL1_chr10.pdf", width = 6, height = 5)
