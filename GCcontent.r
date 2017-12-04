## GENOME GC content
library(data.table)

<<<<<<< HEAD
danRer10_1000<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.danRer.txt")
mm10<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.mm10.txt")
hg38<- fread("/home/emma/mnt/edann/hexamers/GCcont_1000slop_hgTables.refgen.hg38.txt")

danRer10_rand_1000<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.danRer.txt")
mm10.rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.mm10.txt")
hg38.rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand.hg38.txt")


pdf("/home/emma/HexamerBias/output/GC_cont_boxplot.pdf")

boxplot(c(hg38[,4], hg38.rand[,4], mm10[,4],mm10.rand[,4], danRer10[,4], danRer10.rand[,4]),xpd = F, outline = F, names = c("mm10_TSS","mm10_rand", "danRer10_TSS","danRer10_rand", "hg38_TSS", "hg38_rand"), ylab="%GC", main="GC content at TSS (1000 bps flanks)", varwidth = T, col = c("red", "pink"), xaxt="n")
axis(side = 1,at = c(1.5,3.5,5.5),c("hg38", "mm10", "danRer10"),tick = FALSE)
# t.test(hg38[,4], hg38.rand[,4])
# t.test(mm10[,4],mm10.rand[,4])
# t.test(danRer10[,4], danRer10.rand[,4])

legend("topright",c("TSS", "random"), fill=c("red", "pink"), bty="n")

dev.off()

danRer10_400<- fread("/home/emma/mnt/edann/hexamers/GCcont_txStart_danRer10.400.txt")

danRer10_rand_400<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand_danRer10.400.txt")
boxplot(c(danRer10[,4], danRer10_rand[,4]),xpd = F, outline = F, names = c( "danRer10_TSS","danRer10_rand"), ylab="%GC", main="GC content at TSS (400 bps flanks)", varwidth = T, col = c("red", "pink"))

danRer10<- fread("/home/emma/mnt/edann/hexamers/GCcont_txStart_danRer10.200.txt")
danRer10_rand<- fread("/home/emma/mnt/edann/hexamers/GCcont_rand_danRer10.200.txt")

pdf("~/HexamerBias/output/GCcont_danRer.pdf")
boxplot(c(danRer10_1000[,4], danRer10_rand_1000[,4], danRer10_400[,4], danRer10_rand_400[,4], danRer10[,4], danRer10_rand[,4]),xpd = F, outline = F, ylab="%GC", main="GC content at TSS - danRer10", varwidth = T, col = c("red", "pink"), xaxt='n')
legend("topleft",c("TSS", "random"), fill=c("red", "pink"), bty="n")
axis(side = 1,at = c(1.5,3.5,5.5),c("1000 bps", "400 bps", "200 bps"),tick = FALSE)
dev.off()

