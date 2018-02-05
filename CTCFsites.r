## CTCCF sites coverage
library(data.table)

ctcf <- fread("mnt/edann/CTCFdepth_VAN1667.bed")
rand <- fread("mnt/edann/CTCFdepth_rand_VAN1667.bed")
tss <- fread("mnt/edann/TxStart_depth.bed")

pdf("AvOwork/output/CTCFcoverage_wTSS.pdf")
boxplot(tss$V7, c(ctcf$V10, rep(0, 26795 - length(ctcf$V10))), c(rand$V10, rep(0, 26795 - length(rand$V10))), outline = FALSE, varwidth = TRUE, notch = TRUE, names = c('TSS', 'CTCF sites', 'random'), ylab='coverage')
dev.off()
