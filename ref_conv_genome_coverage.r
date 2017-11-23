## Plot coverage datasets
library(data.table)
library(ggplot2)
library(dplyr)

a_chr1<-filter(a,V1=="chr1")
a<-mutate(a,cov = V5+V6)

split3<-lapply(unique(cov3$V1), function(chr) cov3[V1==chr])
split4<-lapply(unique(cov4$V1), function(chr) cov4[V1==chr])

overlap <- sapply(seq_along(split3), function(i) intersect(split3[[i]]$V2,split4[[i]]$V2))

pdf("./output/smp34_cor.pdf", onefile=T)
par(mfrow=c(2,2))
sapply(seq_along(split3), function(i)
	plot(split3[[i]][split3[[i]]$V2 %in% overlap[[i]]]$V4, split4[[i]][split4[[i]]$V2 %in% overlap[[i]]]$V4, 
		pch=19, xlab= "smp 3 met fraction", ylab= "smp 4 met fraction", 
		main=paste(split3[[i]]$V1[1], "( PCC=", cor(split3[[i]][split3[[i]]$V2 %in% overlap[[i]]]$V4, split4[[i]][split4[[i]]$V2 %in% overlap[[i]]]$V4),")") 
		)
	)
dev.off()

sapply(seq_along(split3), function(i)
	cor(split3[[i]][split3[[i]]$V2 %in% overlap[[i]]]$V4, split4[[i]][split4[[i]]$V2 %in% overlap[[i]]]$V4 
		#pch=19, xlab= "smp 3 met fraction", ylab= "smp 4 met fraction", main=split3[[i]]$V1[1] 
		)
	)

## Plot genomic coverage

smp4.cov<-fread("./smp4.genomecov", col.names=c("chr","start","end","cov"), header=F)
plot(smp4.cov[chr=="chr1"]$start, smp4.cov[chr=="chr1"]$cov, type="s")

# split_frac <- lapply(split, function(chr) mutate(chr,frac = V5+V6))

# hist(a$cov, xlim=c(0,200), breaks = 100, main = "CpG site coverage distribution smp.4", xlab="coverage")

## 	HOW MANY UNCOVERED CpGs?
files=list.files("~/mnt/edann/hexamers/kaester/met_extraction/",pattern = "CpG_report", full.names = T)
cov2c_973_1<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454973_1_val_1_bismark_bt2.deduplicated.CpG_report.txt.gz")
cov2c_973_2<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454973_2_val_2_bismark_bt2.deduplicated.CpG_report.txt.gz")
cov2c_965_1<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454965_1_val_1_bismark_bt2.deduplicated.CpG_report.txt.gz")
cov2c_965_2<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454965_2_val_2_bismark_bt2.deduplicated.CpG_report.txt.gz")
cov2c_978_1<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454978_1_val_1_bismark_bt2.deduplicated.CpG_report.txt.gz")
cov2c_978_2<-fread("gunzip -c ~/mnt/edann/hexamers/kaester/met_extraction/ERR454978_2_val_2_bismark_bt2.deduplicated.CpG_report.txt.gz")

cov2c_gen<-fread("gunzip -c ~/mnt/edann/hexamers/cov2c_merged_convRef.txt.gz")

computeCpGcoverage<-function(cov2c){
  cov <- cov2c$V4+cov2c$V5
  tab<-table(cov)
  uncovered<-tab[1]
  covered<-sum(tab[-1])
  coverage<-c(covered,uncovered)
  names(coverage)<-c("covered", "uncovered")
  return(list(coverage,cov))
}

plotCoverage<-function(coverage,cov,smp_name){
  x11()
  par(mfrow=c(1,2))
  pie(coverage,clockwise = T, labels = "", main = smp_name)
  legend("bottomleft",c("covered", "uncovered"), fill=c("lightblue", "white"), bty = "n")
  cov=cov[cov!=0]
  boxplot(cov, ylim=c(0,40), outline = F, main=smp_name, ylab="coverage per CpG")
}
cov2c_beta<- cov2c_gen$V4/(cov2c_gen$V4+cov2c_gen$V5)  


pdf("../output/covered_cpgs_smp3.pdf")
barplot(c(uncovered,covered)/sum(tab), main="smp3 - coverage", names=c("uncovered","covered"), ylim=c(0,1))
dev.off()

cov2c4<-fread("gunzip -c ./hexamers/cov2c_4_tr2_R1_bismark_bt2.deduplicated.bismark.cov.gz")
cov4 <- cov2c4$V4+cov2c4$V5
tab4<-table(cov4)
uncovered<-tab[1]
covered<-sum(tab[-1])

pdf("../output/covered_cpgs_smp1.pdf")
barplot(c(uncovered,covered)/43735674, main="smp4 - coverage", names=c("uncovered","covered"), ylim=c(0,1))
dev.off()

# split_frac <- lapply(split, function(chr) mutate(chr,frac = V5+V6))

# plot(split[[1]]$V2[1:5000], split[[1]]$V4[1:1000]+split[[1]]$V5[1:5000], type="h", border = F)

pdf("./output/hexamers/cov_ref_genomes.pdf")
boxplot(cov,cov3,cov4, ylim=c(0,20), ouline=FALSE, names=c("smp1", "smp3", "smp4"),ylab="coverage per CpG"))
dev.off()