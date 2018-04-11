### AVG METHYLATION: REIK vs NORMAL BS ###
library(data.table)
# source("http://www.bioconductor.org/biocLite.R")
biocLite(c("bsseq"))
library(rtracklayer)
library(devtools)
library(bsseq)
install_github("al2na/methylKit", build_vignettes=FALSE, 
               repos=BiocInstaller::biocinstallRepos(),
               dependencies=TRUE)

load.bismark.met <- function(file, save_file=NULL){
  met <- fread(file, header=FALSE)
  colnames(met) <- c("chr", 'start', 'end', 'frac', 'C', 'T')
  if (!is.null(save_file)) {
    save(met, file = save_file)
  }
  return(met)
  }

reik <- read.bismark("~/mnt/edann/crypts_bs/VAN1667/L1_trim1_R1_bismark_bt2_pe.deduplicated.bismark.thresh5.cov.bed", 'reik', strandCollapse = FALSE)
reik.met <- getMeth(reik, type = 'raw')
wgBS <- read.bismark("~/mnt/edann/hexamers/kaester/met_extraction/ERR454965_1_val_1_bismark_bt2.deduplicated.bismark.thresh5.cov", 'k',strandCollapse = FALSE)
ggplot(reik, aes(frac)) + geom_density(fill='cornflowerblue', alpha=0.2, color='cornflowerblue') +
  theme_classic()
