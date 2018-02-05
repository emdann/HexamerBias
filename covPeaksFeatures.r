### COV PEAKS FEATURES
library(dplyr)
library(ggplot2)
library(data.table)

peaks <- fread('~/multipeaks.nuc.bed')
randpeaks <- fread('mnt/edann/coverage_peaks/multipeaks.random.nuc.bed')

pdf('AvOwork/output/peaks_vs_rand_nuc.pdf')
boxplot(randpeaks$`4_pct_at`, randpeaks$`5_pct_gc`, peaks$`5_pct_at`, peaks$`6_pct_gc`, 
        at = c(1,2,3,4), outline=FALSE, varwidth = TRUE, xaxt='n', ylab='fraction', col = c('cornflowerblue', 'lightpink'))
axis(side = 1, at = c(1.5,3.5), labels=c('random', 'peaks'),tick = FALSE)
legend('topright', c('AT', 'GC'), fill=c('cornflowerblue', 'lightpink'), bty = 'n' )
dev.off()

pdf('AvOwork/output/peaks_vs_rand_CGdinucl.pdf')
boxplot(peaks$`14_user_patt_count`/peaks$`13_seq_len`, randpeaks$`13_user_patt_count`/randpeaks$`12_seq_len`, 
        outline = FALSE, ylab='#CG / length', names = c('peaks', 'random'))
dev.off()
