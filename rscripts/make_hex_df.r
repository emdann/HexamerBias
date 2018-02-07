## MAKE NEW HEX DATAFRAME ##
library(dplyr)
library(caTools)
library(data.table)


hex.df<-read.csv("~/Scaricati/hex_analysis.csv")
usage_df<-readRDS("~/mnt/edann/crypts_bs/VAN1667_hex_usage.RDS")
usage<-apply(usage_df,1,mean)
deltaG<- fread("gunzip -c ~/mnt/edann/hexamers/rand_hex_deltaG.txt.gz")
ab_amp<- read.delim("~/mnt/edann/hexamers/rand_hex_abundance_withamp.txt", header = F, sep= '\t', col.names = c("hex", "abundance_amp"))
ordered_ab_amp<-ab_amp[match(hex.df$hex, ab_amp$hex),]
ab_noamp<- read.delim("~/mnt/edann/hexamers/rand_hex_abundance_noamp.txt", header = F, sep= '\t', col.names = c("hex", "abundance_noamp"))
ordered_ab_noamp<-ab_noamp[match(hex.df$hex, ab_noamp$hex),]

files<-list.files("~/mnt/edann/",pattern = "splitHex", full.names = T)
strand_tabs<-lapply(files, read.delim, header = T, sep = ',')
ord_strand_tabs<-lapply(strand_tabs, function(df) df[match(hex.df$hex, df$X),])
ord_strand_tabs<-lapply(ord_strand_tabs, function(x) mutate(x,OT = ifelse(is.na(OT), 0, OT),OB = ifelse(is.na(OB), 0, OB),CTOT = ifelse(is.na(CTOT), 0, CTOT),CTOB = ifelse(is.na(CTOB), 0, CTOB) ))
strand_df <- data.frame(OB = apply(do.call(cbind,lapply(ord_strand_tabs, function(x) x$OB)), 1,mean),
                        OT = apply(do.call(cbind,lapply(ord_strand_tabs, function(x) x$OT)), 1,mean),
                        CTOB = apply(do.call(cbind,lapply(ord_strand_tabs, function(x) x$CTOB)), 1,mean),
                        CTOT = apply(do.call(cbind,lapply(ord_strand_tabs, function(x) x$CTOT)), 1,mean))
                        
new_df<-data.frame(hex.df$hex, usage=usage, deltaG=deltaG$V2,abundance_noamp=ordered_ab_noamp$abundance_noamp,abundance_amp=ordered_ab_amp$abundance_amp,strand_df)
new_df$X<-NULL
new_df[is.na(new_df)]<-0

write.csv(new_df, file = "~/hex_analysis.csv")
