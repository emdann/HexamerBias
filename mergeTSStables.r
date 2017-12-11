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


TSSdistances<-fread("gunzip -c ~/mnt/edann/hexamers/sumTSS_distances_corrected.txt.gz", header = T)
# TSSdistancesZF<-fread("~/mnt/edann/hexamers/sumTSS_distances_zf.txt",skip = 1, col.names = loadHeaderTSStab("~/mnt/edann/hexamers/TSS_convRef.chr1.txt"))
 
# Make prob function
p_mat <- apply(TSSdistances[,-1],1,function(i) as.numeric(i/sum(i)))
p_df<-cbind(TSSdistances[,1],t(p_mat))
colnames(p_df)<-colnames(TSSdistances)

hex.df<-read.csv("~/hex_analysis.csv")
top_hex<-hex.df[order(hex.df$deltaG, decreasing = F),]$hex[1:7]
bot_hex<-hex.df[order(hex.df$deltaG, decreasing = T),]$hex[1:10]

plotSetProfile<- function(hexs,k=100,title=NULL){
  cols=rainbow(length(hexs))
  idx=which(p_df$hex %in% hexs)
  rm<-apply(p_df[idx,-1],1,runmean,k=k)
  ylim=c(min(rm),max(rm))
  plot(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[idx[1],-1]), k = k), type='l', col=cols[1], ylim=ylim, xlab = "Distance from TSS (bps)", ylab = "Prob")
  sapply(idx[-1], function(hex) lines(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[hex,-1]), k = k), type='l', col=cols[which(idx==hex)]))
  if(!is.null(title)){title(title)}
  }

pdf("~/Van Oudenaarden code/output/TSS_lowest_dG.pdf")
plotSetProfile(top_hex, title = "Lowest Delta G")
plotSetProfile(bot_hex, title = "Lowest Delta G")
dev.off()

# Heatmap
dGsort_hex<-as.character(hex.df[order(hex.df$deltaG),]$hex)
ord_p_df<-p_df[match(dGsort_hex, p_df$hex)]
e<-apply(as.matrix(ord_p_df)[,-1],1, runmean, k=100)
colnames(e)<- ord_p_df$hex

pal=RColorBrewer::brewer.pal(name = "YlOrRd",9)
jpeg("~/Van Oudenaarden code/output/tss_heatmap_mm10.jpg")
heatmap(t(e), Colv = NA, Rowv = NA, labRow =colnames(e), labCol = FALSE, ylab="DeltaG ranked hexamers", revC = TRUE, col=pal, main = "mm10", xlab = "TSS")
# axis(2,outer = T,pos = 0.9, at=seq(-4.5,44,length.out = 5), labels = seq(min(e), max(e), length.out = 5), xpd=T)
# axis(3, at=seq(0.145, 0.725, length.out = 3), labels = c("- 3 kb", "TSS", "3 kb"),outer=TRUE)
dev.off()

Occsort_hex<-as.character(hex.df[order(hex.df$occupancy, decreasing = T),]$hex)
ord_p_df<-p_df[match(Occsort_hex, p_df$hex)]
e<-apply(as.matrix(ord_p_df)[,-1],1, runmean, k=100)
colnames(e)<- ord_p_df$hex

heatmap(t(e), Rowv = NA, Colv = NA, labCol = FALSE, ylab="Occupancy ranked hexamers", col=pal, revC = TRUE)

