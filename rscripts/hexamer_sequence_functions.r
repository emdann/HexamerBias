#####################################################################
###### PRIMER POOL AND HEXAMER SEQUENCE MANIPULATION FUNCTIONS ######
#####################################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(gtools)
# source('~/HexamerBias/rscripts/sanity_checks_deltaG.r')

prevalent_nucleotide <- function(seq){
  nuc.count <- table(strsplit(seq, ''))
  prev.nuc <- names(which.max(nuc.count))
  return(prev.nuc)
}

rev.comp<-function(x,compl=TRUE,rev=TRUE){
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  if(compl==TRUE)
  {
    for (bbb in 1:nchar(x))
    {
      if(xx[bbb]=="A") y[bbb]<-"T"
      if(xx[bbb]=="C") y[bbb]<-"G"
      if(xx[bbb]=="G") y[bbb]<-"C"
      if(xx[bbb]=="T") y[bbb]<-"A"
    }
  }
  if(compl==FALSE)
  {
    for (bbb in 1:nchar(x))
    {
      if(xx[bbb]=="A") y[bbb]<-"A"
      if(xx[bbb]=="C") y[bbb]<-"C"
      if(xx[bbb]=="G") y[bbb]<-"G"
      if(xx[bbb]=="T") y[bbb]<-"T"
    }
  }
  if(rev==FALSE)
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)
}

find.rev.pairs <- function(df,smp){
  p <- filter(df, primer==smp$primer, template==smp$template  | template==rev.comp(smp$template, compl = F) )
  return(p)
}

find.pal.epsilon <- function(smp.ij, smp.ji){
  return(c((smp.ji$abundance*smp.ij$pt - smp.ij$abundance*smp.ji$pt), (smp.ij$pt*smp.ji$cele.pt - smp.ji$pt*smp.ij$cele.pt)))
}

is.palindrome <- function(seq){
  if (seq == paste0(rev(strsplit(seq,split='')[[1]]), collapse='')) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}

## PRIMER POOL SIMULATION ##

build.random.base <- function(pA=0.25, pT=0.25, pG=0.25, pC=0.25){
  a <- pA
  c <- a+pC
  t <- c+pT
  g <- t+pG
  rand <- runif(1,0,1)
  if (rand < a) {return('A')}
  if (rand > a & rand < c) {return('C')}
  if (rand > c & rand < t) {return('T')}
  if (rand > t) {return('G')}
}

build.random.hex <- function(pA=0.25, pT=0.25, pG=0.25, pC=0.25){
  hex <- ''
  for (n in 1:6) {
    hex <- paste0(hex,build.random.base(pA=pA, pT=pT, pG=pG, pC=pC))  
  }
  return(hex)
}

simulate.primer.pool <- function(pool.size=50000, pA=0.25, pT=0.25, pG=0.25, pC=0.25){
  pool <- c()
  for (n in 1:pool.size) {
    pool <- c(pool, build.random.hex(pA=pA, pT=pT, pG=pG, pC=pC))
  }
  return(pool)
}

all.hexamers <- function(){
  hexs <- apply(permutations(4,6, v=c("A", "C", "T", "G"), repeats.allowed = T),1, function(x) paste(x, collapse = ""))
  return(hexs)
  }

primer.prob <- function(seq, probs = c(pA=0.25, pT=0.25, pG=0.25, pC=0.25)){
  prob=1
  for (l in strsplit(seq, '')) {
    prob <- prob*probs[paste0("p",l)]
  }
  return(prod(prob))
}

batch.prob.uniform <- function(hexs = all.hexamers(), nuc.probs = c(pA=0.25, pT=0.25, pG=0.25, pC=0.25)){
  # Computes probability for each hexamer given a set of probabilities for each nucleotide
  # N.B. SAME PROBABILITY FOR ALL POSITIONS!
  probs <- sapply(hexs, primer.prob, probs=nuc.probs)
  return(probs)
  }

hexamerMatrix <- function(stepSize = 0.1){
  # Computes all possible probability combinations for 4 nucletides given a step size for the difference in probability
  vals = seq(from = 0, to = 1, by = stepSize)
  combMat = expand.grid(vals, vals, vals)
  rowSums = apply(combMat, 1, sum)
  impVals = rowSums > 1
  combMat = combMat[!impVals, ]
  combMat[, 4] = (1 - apply(combMat, 1, sum))
  colnames(combMat) = c("pA", "pC", "pT", "pG")
  return(combMat)
}

