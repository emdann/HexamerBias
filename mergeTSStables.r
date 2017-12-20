library(dplyr)
library(data.table)
library(caTools)
library(RColorBrewer)
library(gplots)
library(heatmap3)
library(pheatmap)


TESdistances<-fread("~/mnt/edann/hexamers/sumTES_distances_conv.txt", header = T)
TSSdistances<-fread("~/mnt/edann/hexamers/sumTSS_distances_recorrected.txt", header = T)
Rand_distances<-fread("~/mnt/edann/hexamers/sumRand_conv_test.txt", header = T)
CTCFdistances<-fread("~/mnt/edann/hexamers/annotations_bed/CTCFcount_conv.txt", header = T)
Exstart_distances<-fread("~/mnt/edann/hexamers/annotations_bed/exonStarts_count_conv.txt", header = T)

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

makeDeltaGHeatmap <- function(prob, center="TSS"){
  dGsort_hex<-as.character(hex.df[order(hex.df$deltaG),]$hex)
  ord_p_df<-prob[match(dGsort_hex, prob$hex)]
  ord_p_nona <- ord_p_df[!apply(ord_p_df,1, function(x) all(is.na(x))),]
  e<-apply(as.matrix(ord_p_nona)[,-1],1, runmean, k=100)
  colnames(e)<- ord_p_nona$hex
  pal=colorRampPalette(brewer.pal(name = "YlOrRd",9))(5000)
  s <- seq(-2999, 2995) 
  anno_row <- data.frame(`deltaG(kcal/mol)`=sort(hex.df$deltaG), row.names = hex.df[order(hex.df$deltaG),]$hex)
  pheatmap(t(e), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, annotation_row=anno_row, show_rownames = FALSE, annotation_legend = TRUE , show_colnames = TRUE,labels_col = ifelse(s %in% c(-2999, 0, 2995),s,''), col=pal, fontsize_row = 1, filename = "~/Van Oudenaarden code/output/pretty_heatmap.png")
  # heatmap(t(e), col=pal, xlab = center)
}

tes_p <- compute_prob(TESdistances)
tss_p <- compute_prob(TSSdistances)
rand_p <- compute_prob(Rand_distances)
ctcf_p <- compute_prob(CTCFdistances)
exstart_p<-compute_prob(Exstart_distances)


hex.df<-read.csv("~/hex_analysis.csv")
top_hex<-hex.df[order(hex.df$deltaG, decreasing = F),]$hex[1:10]
bot_hex<-hex.df[order(hex.df$deltaG, decreasing = T),]$hex[1:10]

pdf("~/Van Oudenaarden code/output/TES_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(tes_p, top_hex, center="TES", title = "Lowest Delta G", leg.pos = "topright")
plotSetProfile(tes_p, bot_hex, center="TES", title = "Highest Delta G")
dev.off()

pdf("~/Van Oudenaarden code/output/Rand_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(rand_p,top_hex, center= "random site", title = "Lowest Delta G")
plotSetProfile(rand_p,bot_hex, center= "random site", title = "Highest Delta G")
dev.off()

pdf("~/Van Oudenaarden code/output/Rand_topbot_dG.pdf", width=12)
par(mfrow=c(1,2))
plotSetProfile(tes_p,top_hex, center= "TES", title = "Lowest Delta G", leg.pos = "topright")
plotSetProfile(tes_p,bot_hex, center= "TES", title = "Highest Delta G")
dev.off()



# Heatmap
makeDeltaGHeatmap(rand_p)
x11()
makeDeltaGHeatmap(tes_p)
x11()
makeDeltaGHeatmap(tss_p)


# heatmap3(t(e),Rowv = NA,Colv = NA, showColDendro = FALSE,showRowDendro = FALSE, labRow =FALSE, labCol = FALSE, ylab="DeltaG ranked hexamers", revC = TRUE, col=pal, main = "mm10", xlab = "TSS")
# axis(2,outer = T,pos = 0.9, at=seq(-4.5,44,length.out = 5), labels = seq(min(e), max(e), length.out = 5), xpd=T)
# axis(3, at=seq(0.145, 0.725, length.out = 3), labels = c("- 3 kb", "TSS", "3 kb"),outer=TRUE)
dev.off()

png("~/Van Oudenaarden code/output/tss_heatmap_unconv_mm10.png", width = 550, height = 550)
Occsort_hex<-as.character(hex.df[order(hex.df$occupancy, decreasing = T),]$hex)
ord_p_df<-p_df[match(Occsort_hex, p_df$hex)]
e<-apply(as.matrix(ord_p_df)[,-1],1, runmean, k=100)
colnames(e)<- ord_p_df$hex

heatmap(t(e), Rowv = NA, Colv = NA, labCol = FALSE, ylab="Occupancy ranked hexamers", col=pal, revC = TRUE)

