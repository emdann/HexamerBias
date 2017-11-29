library(dplyr)
library(data.table)
library(caTools)
library(RColorBrewer)


loadHeaderTSStab<-function(file_path){
  header<-scan(file_path,what = "", nlines = 1)
  hex<-substr(header[1],1,4)
  header[1]<-substr(header[1],5,11)
  return(header)
}

TSSdistances<-fread("~/mnt/edann/hexamers/sumTSS_distances.txt", col.names = c("hex", loadHeaderTSStab("~/mnt/edann/hexamers/TSS_tab.chr10.txt")))

# Make prob function
p_mat <- apply(t(TSSdistances[,-1]),1,function(i) as.numeric(i/sum(i)))
p_df<-cbind(TSSdistances[,1],p_mat)

sapply(sample(seq_along(p_df$hex),10), function(hex) lines(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[hex,-1]), k = 100), type='l', main=paste(p_df$hex[hex])))


hex.df<-read.csv("~/hex_analysis.csv")
top_hex<-hex.df[order(hex.df$deltaG, decreasing = F),]$hex[1:10]
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

pdf("~/HexamerBias/output/TSS_profiles_mostVSleastused.pdf", width = 12)
par(mfrow=c(1,2))
plotSetProfile(top_hex, title = "Most used hexamers")
plotSetProfile(bot_hex, title = "Least used hexamers")
dev.off()

hclustfunc <- function(x, method = "complete", dmeth = "euclidean") {    
  hclust(dist(x, method = dmeth), method = method)
}

