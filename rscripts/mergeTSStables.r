library(dplyr)
library(data.table)
library(caTools)
library(RColorBrewer)
library(gplots)
library(heatmap3)
library(pheatmap)
library(rPython)
python.load("/home/emma/HexamerBias/freeEnergy.py")

# Make prob function
compute_prob <- function(dist_matrix){
  p_mat <- apply(dist_matrix[,-1],1,function(i) as.numeric(i/sum(i)))
  p_df<-cbind(dist_matrix[,1],t(p_mat))
  colnames(p_df)<-colnames(dist_matrix)
  return(p_df)
}

plotSetProfile<- function(p_df,hexs,k=100,title=NULL, center="TSS", leg.pos="bottomright"){
  cols=rainbow(length(hexs))
  idx=which(p_df$hex %in% hexs)
  rm<-apply(p_df[idx,-1],1,runmean,k=k)
  ylim=c(min(rm),max(rm))
  plot(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[idx[1],-1]), k = k), type='l', col=cols[1], ylim=ylim, xlab = paste("Distance from",center, "(bps)"), ylab = "Prob")
  sapply(idx[-1], function(hex) lines(as.numeric(colnames(p_df)[-1]),runmean(as.numeric(p_df[hex,-1]), k = k), type='l', col=cols[which(idx==hex)]))
  if(!is.null(title)){title(title)}
  legend(leg.pos, legend = hexs, lty = 1, col=cols, bty='n', cex=.6)
}

# makeDeltaGHeatmap <- function(prob, center="TSS"){
#   dGsort_hex<-as.character(hex.df[order(hex.df$deltaG),]$hex)
#   ord_p_df<-prob[match(dGsort_hex, prob$hex)]
#   ord_p_nona <- ord_p_df[!apply(ord_p_df,1, function(x) all(is.na(x))),]
#   e<-apply(as.matrix(ord_p_nona)[,-1],1, runmean, k=100)
#   colnames(e)<- ord_p_nona$hex
#   pal=colorRampPalette(brewer.pal(name = "YlOrRd",9))(5000)
#   s <- seq(-2999, 2995) 
#   anno_row <- data.frame(`deltaG(kcal/mol)`=sort(hex.df$deltaG), row.names = hex.df[order(hex.df$deltaG),]$hex)
#   pheatmap(t(e), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, annotation_row=anno_row, show_rownames = FALSE, annotation_legend = TRUE , show_colnames = TRUE,labels_col = ifelse(s %in% c(-2999, 0, 2995),s,''), col=pal, fontsize_row = 1, filename = "~/Van Oudenaarden code/output/pretty_heatmap.png")
#   # heatmap(t(e), col=pal, xlab = center)
# }

# Different K-lengths
makeDeltaGHeatmap <- function(tss_p, k){
  noN_p <- tss_p[-grep("N", tss_p$hex)]
  sort_dG<- sort(sapply(noN_p$hex, function(seq) python.call("compute_deltaG_ions", seq)), decreasing = F)
  e<-apply(as.matrix(noN_p[match(names(sort_dG), noN_p$hex)])[,-1],1, runmean, k=100)
  colnames(e)<- names(sort_dG)
  pal=colorRampPalette(brewer.pal(name = "YlOrRd",9))(5000)
  s <- seq(-2999, 2995) 
  anno_row <- data.frame(`deltaG(kcal/mol)`=sort_dG, row.names = colnames(e))
  pheatmap(t(e), cluster_rows = FALSE, cluster_cols = FALSE , border_color = NA, 
           annotation_row=anno_row, show_rownames = FALSE, annotation_legend = TRUE, 
           show_colnames = TRUE, labels_col = ifelse(s %in% c(-2999, 0, 2995),s,''), 
           col=pal, main=paste("k = ",k), filename = paste0("~/Van Oudenaarden code/output/heatmap_k",k,".pdf"))
}


TSSdistances<-fread("~/mnt/edann/hexamers/sumTSS_k6_unconv.txt", header = T, na.strings = ",,")
tss_p <- compute_prob(TSSdistances)
makeDeltaGHeatmap(tss_p, 6)

hex.df<-read.csv("~/hex_analysis.csv")
top_hex<-hex.df[order(hex.df$deltaG, decreasing = F),]$hex[1:10]
bot_hex<-hex.df[order(hex.df$deltaG, decreasing = T),]$hex[1:10]

pdf("~/Van Oudenaarden code/output/ExStart_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(exstart_p, top_hex, center="Exon start", title = "Lowest Delta G", leg.pos = "topright")
plotSetProfile(exstart_p, bot_hex, center="Exon start", title = "Highest Delta G")
dev.off()

pdf("~/Van Oudenaarden code/output/Rand_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(tss_p,top_hex, center= "TSS", title = "Lowest Delta G")
plotSetProfile(tss_p,bot_hex, center= "TSS", title = "Highest Delta G")
dev.off()

pdf("~/Van Oudenaarden code/output/Rand_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(tes_p,top_hex, center= "TES", title = "Lowest Delta G", leg.pos = "topright")
plotSetProfile(tes_p,bot_hex, center= "TES", title = "Highest Delta G")
dev.off()



# Heatmap



