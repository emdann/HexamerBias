## Plot coverage datasets
library(data.table)
library(ggplot2)
library(RColorBrewer)

colpal=RColorBrewer::brewer.pal(10,"Set3")

## 	HOW MANY UNCOVERED CpGs?
files=list.files("/hpc/hub_oudenaarden/edann/hexamers/kaester/met_extraction/",pattern = "CpG_report", full.names = T)

cov2c_gen<-fread("gunzip -c ~/mnt/edann/hexamers/cov2c_merged_convRef.txt.gz")

computeCpGcoverage<-function(f){
  cov2c<-fread(paste("gunzip -c", f))
  cov <- cov2c$V4+cov2c$V5
  tab<-table(cov)
  uncovered<-tab[1]
  covered<-sum(tab[-1])
  coverage<-c(covered,uncovered)
  names(coverage)<-c("covered", "uncovered")
  cov2c_beta<- cov2c$V4/(cov2c$V4+cov2c$V5) 
  gc()
  return(list(coverage,cov,cov2c_beta))
}


cov2c<-sapply(files, computeCpGcoverage)
smp_names=sapply(files, function(f) gsub(pattern = "_val.+", replacement = "",strsplit(f,"//")[[1]][2]))
ref<-computeCpGcoverage("/hpc/hub_oudenaarden/edann/hexamers/cov2c_merged_convRef.txt.gz")

plotCoverage<-function(coverage,cov,smp_name){
  x11()
  par(mfrow=c(1,2))
  pie(coverage,clockwise = T, labels = "", main = smp_name)
  legend("bottomleft",c("covered", "uncovered"), fill=c("lightblue", "white"), bty = "n")
  cov=cov[cov!=0]
  boxplot(cov_dist[[1]],cov_dist[[2]],cov_dist[[3]],cov_dist[[4]],cov_dist[[5]],cov_dist[[6]],cov_dist[[7]], ylim=c(0,40), outline = F, main=smp_name, ylab="coverage per CpG")
}
 
cov_dist=append(cov2c[c(2,5,8,11,14,17)], ref)[1:8][-7]
pdf("/hpc/hub_oudenaarden/edann/output/hexamers/refConvGenome_stats/cov_dist.pdf")
par(mar=c(7,4,2,2)+0.2)
boxplot(cov_dist[[1]],cov_dist[[2]],cov_dist[[3]],cov_dist[[4]],cov_dist[[5]],cov_dist[[6]],cov_dist[[7]], 
        ylim=c(0,40), outline = F, ylab="coverage per CpG", names = c(smp_names, "averaged"), las=2)
dev.off()
