## GENOME GC content
library(data.table)

danRer10<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.danRer.txt")
mm10<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.mm10.txt")
hg38<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.hg38.txt")

danRer10.rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.danRer.txt")
mm10.rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.mm10.txt")
hg38.rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.hg38.txt")


pdf("/home/emma/HexamerBias/output/GC_cont_boxplot.pdf")

boxplot(c(hg38[,4], hg38.rand[,4], mm10[,4],mm10.rand[,4], danRer10[,4], danRer10.rand[,4]),xpd = F, outline = F, names = c("mm10_TSS","mm10_rand", "danRer10_TSS","danRer10_rand", "hg38_TSS", "hg38_rand"), ylab="%GC", main="GC content at TSS (1000 bps flanks)", varwidth = T, col = c("red", "pink"), xaxt="n")
axis(side = 1,at = c(1.5,3.5,5.5),c("hg38", "mm10", "danRer10"),tick = FALSE)
# t.test(hg38[,4], hg38.rand[,4])
# t.test(mm10[,4],mm10.rand[,4])
# t.test(danRer10[,4], danRer10.rand[,4])
legend("topright",c("TSS", "random"), fill=c("red", "pink"))

dev.off()
