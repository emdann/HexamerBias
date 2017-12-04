library(dplyr)
library(data.table)
library(caTools)
library(RColorBrewer)
library(gplots)



loadHeaderTSStab<-function(file_path){
  header<-scan(file_path,what = "", nlines = 1)
  hex<-substr(header[1],1,4)
  header[1]<-substr(header[1],5,11)
  return(header)
}


TSSdistances<-fread("gunzip -c ~/mnt/edann/hexamers/sumTSS_distances.txt.gz", col.names = c("hex",seq(-2995,2999)))

# Make prob function
p_mat <- apply(t(TSSdistances[,-1]),1,function(i) as.numeric(i/sum(i)))
p_df<-cbind(TSSdistances[,1],p_mat)

hex.df<-read.csv("~/hex_analysis.csv")
top_hex<-hex.df[order(hex.df$deltaG, decreasing = F),]$hex[1:10]
bot_hex<-hex.df[order(hex.df$deltaG, decreasing = T),]$hex[1:10]


top_ab_hex<-hex.df[order(hex.df$abundance, decreasing = T),]$hex[1:10]
bot_ab_hex<-hex.df[order(hex.df$abundance, decreasing = F),]$hex[1:10]

pdf("~/HexamerBias/output/TSS_profiles_mostVSleastabund.pdf", width = 12)
par(mfrow=c(1,2))
plotSetProfile(c(as.character(top_hex),as.character(bot_hex)), title = "Usage")
plotSetProfile(bot_ab_hex, title = "Least abundant hexamers")
dev.off()


plotSetProfile<- function(hexs,k=100,title=NULL){
  cols=rainbow(length(hexs))
  idx=which(p_df$hex %in% hexs)
  rm<-apply(p_df[idx,-1],1,runmean,k=k)
  ylim=c(min(rm),max(rm))
  plot(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[idx[1],-1]), k = k), type='l', col=cols[1], ylim=ylim, xlab = "Distance from TSS (bps)", ylab = "Prob")
  sapply(idx[-1], function(hex) lines(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[hex,-1]), k = k), type='l', col=cols[which(idx==hex)]))
  if(!is.null(title)){title(title)}
  }

pdf("~/HexamerBias/output/TSS_profiles_mostVSleastused.pdf", width = 12)
par(mfrow=c(1,2))
plotSetProfile(top_hex, title = "Most used hexamers")
plotSetProfile(bot_hex, title = "Least used hexamers")
dev.off()


# Heatmap
dGsort_hex<-as.character(hex.df[order(hex.df$deltaG),]$hex)
ord_p_df<-p_df[match(dGsort_hex, p_df$hex)]
e<-apply(as.matrix(ord_p_df)[,-1],1, runmean, k=100)
colnames(e)<- ord_p_df$hex

pdf("~/HexamerBias/output/tss_heatmap_mm10.pdf")
heatmap(t(e), Rowv = NA, Colv = NA, labRow =FALSE, labCol = FALSE, ylab="DeltaG ranked hexamers", col=pal,revC = TRUE, RowSideColors = pal2)
# axis(2,outer = T,pos = 0.9, at=seq(-4.5,44,length.out = 5), labels = seq(min(e), max(e), length.out = 5), xpd=T)
axis(1, pos=-0.1, at=seq(0.14, 0.77, length.out = 3), labels = c("- 3 kb", "TSS", "3 kb"))
dev.off()

Occsort_hex<-as.character(hex.df[order(hex.df$occupancy, decreasing = T),]$hex)
ord_p_df<-p_df[match(Occsort_hex, p_df$hex)]
e<-apply(as.matrix(ord_p_df)[,-1],1, runmean, k=100)
colnames(e)<- ord_p_df$hex

heatmap(t(e), Rowv = NA, Colv = NA, labCol = FALSE, ylab="Occupancy ranked hexamers", col=pal, revC = TRUE)

