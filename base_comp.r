## base composition of reads
# How did I make the files? 
# zcat CpG_OT_L2_trim1_R1_bismark_bt2_pe.deduplicated.txt.gz | head -1000 |  cut -f 1 > OT_IDS.test
# /hpc/hub_oudenaarden/bdebarbanson/bin/samtools-1.4.1/samtools fastq L2_trim1_R1_*deduplicated.bam | grep -f OT_IDS.test -A1 | grep -v "^@" | grep -v - > seqs_OT_L2.test

files_or<-list.files("~/mnt/organoids_bs/split_orgs_1", pattern = "original_.+.hextab", full.names = T)
original_fastas_spl1<-lapply(files_or, readLines)
files_compl<-list.files("~/mnt/organoids_bs/split_orgs_1", pattern = "compl_.+.hextab", full.names = T)
compl_fastas_spl1<-lapply(files_compl, scan, what="")


base_comp<- function(s){
  countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))}
  t <- countCharOccurrences("T",s)/nchar(s)
  a<-countCharOccurrences("A",s)/nchar(s)
  g<-countCharOccurrences("G",s)/nchar(s)
  c<-countCharOccurrences("C",s)/nchar(s)
  return(c(a,t,c,g))
}
base_comp_CTOB <- do.call(rbind, tapply(CTOB, seq_along(CTOB), base_comp))
barplot(t(base_comp_CTOB)[,1:50], col = rainbow(4), legend.text = rownames(t(base_comp_CTOB)))

base_comp_OT <- do.call(rbind, tapply(OT, seq_along(OT), base_comp))
barplot(t(base_comp_OT)[,1:50], col = rainbow(4), legend.text = rownames(t(base_comp_OT)))


plot(cumsum(sort(table(substr(compl,1,6)), decreasing=T)), pch=".")
