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


## Plot percentage of mismatcing primer-template events 

load("~/AvOwork/rdata/cele_pt_all.RData")
mm.count <- cele.all.df %>% 
  mutate(sample='cele_noBS') %>%
  bind_rows(., mutate(human.all.df, sample='human_noBS')) %>%
  bind_rows(., mutate(zf.all.df, sample='zf_noBS')) %>%
  bind_rows(., mutate(d3r.all.df, sample='D3R_BS')) %>%
  mutate(MM=ifelse(template==primer, 'match', 'mismatch')) %>% 
  group_by(sample, MM) %>%
  summarise(pt=sum(pt)) 
  
mm.count %>%
  filter(sample!="D3R_BS") %>%
  ggplot(., aes(sample, pt*100, fill=MM)) + 
  theme_classic() +
  geom_bar(position='fill', stat='identity') +
  ylab("% Reads") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=30)) +
  scale_fill_brewer(palette = "Reds")

ggsave("~/AvOwork/output/mismatch_noBS.pdf")
