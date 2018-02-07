library(ggplot2)
library(dplyr)
library(reshape2)

L1<- read.delim("~/mnt/edann/L2.splitHex", header = T, sep = ',')
L1[is.na(L1)]<-0

p<- L1 %>% mutate(tot = OT+CTOT+CTOB+OB) %>% filter(., !grepl("N", X)) %>% 
  mutate(ratioOC=(OT+OB)/(CTOB+CTOT)) %>% 
  arrange(.,ratioOC) #%>%
#   ggplot(., aes(x=X, ratioOC, ordered=FALSE)) +
#   geom_bar(stat="identity") + xlab(label = NULL)
# # #  scale_x_continuous(breaks = .$hex, labels = NULL)

p <- p %>% mutate(O=OT+OB, C=CTOB+CTOT)

hex.df<-read.csv("~/hex_analysis.csv")
rand_hex_abundance<- read.delim("~/mnt/edann/hexamers/rand_hex_abundance_noamp.txt", header = F, sep= '\t', col.names = c("hex", "abundance"))
ordered_ab_noamp<-rand_hex_abundance[match(hex.df$hex, rand_hex_abundance$hex),]

pdf("~/Van Oudenaarden code/output/occupancy_deltaG_OriginalUsage.pdf")
plot(hex.df$deltaG,log(p$O/ordered_ab_noamp$abundance), pch=19, ylab="log Occupation no. (usage/abundance)", xlab="DeltaG", main = "Usage on original strand")
text(hex.df$deltaG, log(p$O/ordered_ab_noamp$abundance)+0.2, hex.df$hex, cex=0.6)
dev.off()
