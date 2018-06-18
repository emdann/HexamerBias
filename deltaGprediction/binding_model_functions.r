### Keq COMPUTATION ###
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(fitdistrplus)
source('~/HexamerBias/rscripts/sanity_checks_deltaG.r')

# Putting together functions used in normalization_strategy.Rmd



## DATA PARSING ##

loadPtMatrix <- function(file, compression='gzip'){
  if (compression=='gzip') {
    tab <- fread(paste('zcat',file), sep = ',', header=TRUE)
  }else{
    tab <- fread(file, sep = ',', header=TRUE) 
  }
  p <- as.matrix(tab[,-1])
  rownames(p)<- tab$V1
  return(p)
}

make.hex.usage.df <- function(ptTab, type = 'primer', scale =TRUE){
  if (type=='primer') { hex.usage <- colSums(ptTab)  }
  if (type=='template') { hex.usage <- rowSums(ptTab)  }
  sample.name <- gsub(deparse(substitute(ptTab)), pattern = 'pt.', replacement = '')
  if (scale) { 
    hex.df <- data.frame(names(hex.usage), hex.usage/sum(hex.usage), row.names = NULL) 
  }else{
    hex.df <- data.frame(names(hex.usage), hex.usage, row.names = NULL)
  }
  colnames(hex.df) <- c(type,sample.name)
  return(hex.df)
}

load.pt.data <- function(ptCounts.file, diag.pairs = F){
  # Loading primer-template matrix, returns long matrix 
  pt <- loadPtMatrix(ptCounts.file,compression = 'none')
  print("primer-template matrix loaded!")
  template.usage <- make.hex.usage.df(pt, type='template', scale=F)
  print("template usage dataframe done!")
  primer.usage <- make.hex.usage.df(pt, type='primer', scale=F)
  print("primer usage dataframe done!")
  colnames(template.usage) <- c('template', 't.usage')
  colnames(primer.usage) <- c('primer', 'p.usage')
  pairs <- make_pair_df(pt)
  print("Long dataframe for pairs done!")
  if (diag.pairs) {
    matches <- pairs[substr(pairs$ptPair,1,6)==substr(pairs$ptPair,8,13),]
  } else {
    matches <- pairs
  }
  return(list(matches=matches, t.usage=template.usage, p.usage=primer.usage))
}

load.kmer.abundance <- function(kmer.abundance.file){
  return(read.csv(kmer.abundance.file, header = F, col.names = c('template', 'abundance')))
}

load.modelled.deltaG <- function(deltaG.file){
  return(read.csv(deltaG.file, header = F, col.names = c('template', 'dG'), sep=' '))
}

join.pt.data <- function(matches, template.usage, kmer.ab, tabDg){
  df <- matches %>%
    # filter(dG!=0) %>%
    transmute(template=substr(matches$ptPair,1,6), primer=substr(matches$ptPair,8,100), pt=dG) %>%
    inner_join(., kmer.ab, by='template') %>%
    inner_join(., template.usage, by='template') %>%
    inner_join(., tabDg, by='template')
  return(df)
}

select.diag.pairs <- function(pairs.df){
  pairs.df <- filter(pairs.df, template == primer)
  return(pairs.df)
}

make.match.df <- function(pt.file, ab.file){
  pt.data <- load.pt.data(ptCounts.file = pt.file, diag.pairs = F)
  kmer.ab.data <- load.kmer.abundance(ab.file)
  deltaG.data <- load.modelled.deltaG("~/mnt/edann/hexamers/rand_hex_deltaG_ions.txt.gz")
  matches <- filter(pt.data$matches, substr(ptPair,1,6)==substr(ptPair, 8,100))
  hex.df <- join.pt.data(matches, pt.data$t.usage, kmer.ab.data, deltaG.data)
  return(hex.df)
}

## NORMALIZATION ## 
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

epsilon.distribution <- function(hex.df){
  pals <- unique(hex.df$primer)[sapply(unique(hex.df$primer), is.palindrome)]
  pal.df <- hex.df %>%
    filter(primer %in% pals)
  l.pairs.pals <- lapply(1:nrow(pal.df), function(i) find.rev.pairs(pal.df,pal.df[i,]))
  epsilons.pals <- sapply(l.pairs.pals, function(p) find.pal.epsilon(p[1,], p[2,]))
  return(epsilons.pals)
}

# eps <- median(epsilons, na.rm=T)

normalize.w.palindrome.primers <- function(hex.df, eps, C=1){
  scaled.hex.df <- hex.df %>%
    mutate(scaled.pt = pt*eps*C) %>%
    mutate(pred.dg.scaled= (abundance-scaled.pt)/scaled.pt) %>%
    mutate(pred.dg= (abundance-pt)/pt) 
  return(scaled.hex.df)
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

build.random.hex <- function(){
  hex <- ''
  for (n in 1:6) {
    hex <- paste0(hex,build.random.base())  
  }
  return(hex)
}

simulate.primer.pool <- function(pool.size=50000){
  pool <- c()
  for (n in 1:pool.size) {
    pool <- c(pool, build.random.hex())
  }
  return(pool)
}

# tab <- data.frame(pool = factor())
# for (n in 1:100) {
#   pool <- simulate.primer.pool(10000)
#   count <- as.data.frame(table(pool))
#   tab <- full_join(tab, count, by='pool')
# }

## SCALING FACTOR ESTIMATION ASSUMING DISTRIBUTION OF PRIMER CONCENTRATION ##
group.coeffs.wPrimer.all <- function(hex.df, groupsOI, gamma.shape=24.05, gamma.rate=0.98, pool.size=sum(table(pool1)), imposeP=F){
  groups.df <- hex.df %>% 
    inner_join(., t.groups, by='template') %>%
    filter(group %in% groupsOI) %>%
    mutate(dG=dG/0.59) 
  if (imposeP) {
    groups.df <- mutate(groups.df, p=1/(4^6))
  } else {
    groups.df <- mutate(groups.df, 
                        p=(sample(rgamma(pool.size, 
                                         shape=gamma.shape, 
                                         rate=gamma.rate)/pool.size, n())))
  }
  groups.df.coeffs <- groups.df %>%
    mutate(frac.abundance=(abundance/sum(abundance))) %>%
    mutate(A=sum((p^2)*((frac.abundance/pt)^2)),
           C=sum(((frac.abundance*t.usage)/(pt^2))*(p^2))) %>%
    group_by(group) %>%
    mutate(B=sum((p*frac.abundance)/pt),
           D=sum((t.usage/pt)*p),
           N=n(),
           C.sing=sum(((frac.abundance*t.usage)/(pt^2))*(p^2))
    ) %>%
    dplyr::select(A,B,D,N,C,dG) %>%
    summarise_all(funs(first)) 
  return(list(df=groups.df, coeffs=groups.df.coeffs))  
}

solve.system <- function(groups.df){
  first.left <- c(groups.df$A[1], groups.df$B)
  first.right <- groups.df %>%
    transmute(C) %>%
    sample_n(1)
  first.right <- c(eps=first.right$C)
  coeffs <-c()
  for (n in groups.df$group){
    coeffs <- rbind(coeffs, sapply(groups.df$group, function(i) ifelse(groups.df[groups.df$group==i,]$group==n,as.numeric(groups.df[groups.df$group==i,'N']), 0)))
    # ifelse(groups.df$group==n)
    }
  right <- groups.df %>%
    transmute(group,D)
  right <- as.matrix(cbind(right$group, right$D))
  rownames(right) <- right[,1]
  right <- right[,2]
  left.mat <- rbind(as.numeric(first.left), cbind(groups.df$B, coeffs))
  right.mat <- c(first.right, right)
  # print(left.mat)
  if (any(is.infinite(left.mat[,1]))) {
    return(NA)
  } else {
    # left.mat <- left.mat[!is.infinite(left.mat[,1]),]
    # print(left.mat)
    # print(Det(left.mat))
    sols <- solve(left.mat, right.mat, tol = 1e-18)
    eps <- sols[1]
    dgs <- sols[2:length(sols)]
    return(list(eps=eps, dgs=dgs))  
  }
}

estimate.epsilon <- function(hex.df, primer.pool, groupsOI, imposeP = F, iterations=1000){
  fit.gamma <- summary(fitdist(as.numeric(table(primer.pool)), 'gamma'))
  epsilons <- c()
  i <- 0
  while(i<iterations){
    i<- i+1
    coeffs <- group.coeffs.wPrimer.all(hex.df, groupsOI = groupsOI, gamma.shape = fit.gamma$estimate['shape'],
                                       gamma.rate = fit.gamma$estimate['rate'],
                                       pool.size = length(primer.pool),
                                       imposeP=imposeP)
    sols <- solve.system(coeffs$coeffs)
    if(!is.na(sols)){
      epsilons <- c(epsilons, sols$eps)
    } else {epsilons <- c(epsilons, NA)}
    }
  return(epsilons)
  }

estimate.keq <- function(hex.df, primer.pool, groupsOI, iterations=1000){
  fit.gamma <- summary(fitdist(as.numeric(table(primer.pool)), 'gamma'))
  keq <- c()
  i <- 0
  while(i<iterations){
    i<- i+1
    coeffs <- group.coeffs.wPrimer.all(hex.df, groupsOI = groupsOI, gamma.shape = fit.gamma$estimate['shape'],
                                       gamma.rate = fit.gamma$estimate['rate'],
                                       pool.size = length(primer.pool))
    sols <- solve.system(coeffs$coeffs)
    if(!is.na(sols)){
      keq <- rbind(keq, sols$dgs)
    } else {keq <- rbind(keq, rep(NA, length(groupsOI)))}
  }
  keq.df <- data.frame(keq)
  colnames(keq.df) <- as.factor(groupsOI)
  return(keq.df)
}

estimate.deltaG <- function(hex.df, primer.pool, groupsOI, iterations=1000){
  fit.gamma <- summary(fitdist(as.numeric(table(primer.pool)), 'gamma'))
  deltaG <- c()
  i <- 0
  while(i<iterations){
    i<- i+1
    coeffs <- group.coeffs.wPrimer.all(hex.df, groupsOI = groupsOI, gamma.shape = fit.gamma$estimate['shape'],
                                       gamma.rate = fit.gamma$estimate['rate'],
                                       pool.size = length(primer.pool))
    sols <- solve.system(coeffs$coeffs)
    keq <- sols$dgs
    dgs.df <- as.data.frame(keq) %>%
      mutate(group=groupsOI) %>%
      # rename(Keq=sols$dg) %>%
      inner_join(., coeffs$df) %>%
      mutate(deltaG=keq/p)
  #   if(!is.na(sols)){
  #     keq <- rbind(keq, sols$dgs)
  #   } else {keq <- rbind(keq, rep(NA, length(groupsOI)))}
    deltaG <- rbind(deltaG, dgs.df$deltaG)
    }
  deltaG.df <- data.frame(deltaG)
  colnames(deltaG.df) <- dgs.df$template
  return(deltaG.df)
}

## DeltaG computation

# function(hex.df, primer.pool, epsilon)
# 


 
# fit.gamma <- summary(fitdist(as.numeric(table(primer.pool)), 'gamma'))
# gamma.shape = fit.gamma$estimate['shape']
# gamma.rate = fit.gamma$estimate['rate']
# pool.size = length(primer.pool)
# 
add.many.ps <- function(df, gamma.shape, gamma.rate, pool.size, n.iterations=10){
  start.cols <- colnames(df)
  i <- 0
  while (i<n.iterations){
    i<- i+1
    df <- mutate(df, p=(sample(rgamma(pool.size,shape=gamma.shape,rate=gamma.rate)/pool.size, n())))
    colnames(df)[ncol(df)] <- paste0('p',i)
  }
  df.long <- df %>%
    melt(id.vars=start.cols, variable.name='p.iter', value.name = 'p')
  return(df.long)
}

compute.keqs <- function(pt.dfs = list(cele=cele.df, human=human.df, zf=zf.df), cele.eps, human.eps, zf.eps, gamma.shape, gamma.rate, pool.size, n.iterations=10, take.pairs=F){
  if (take.pairs) {
    dfs.w.p <- lapply(pt.dfs, add.pairs.info)
  } else {
    dfs.w.p <- pt.dfs
  }
  dfs.w.p <- lapply(dfs.w.p, function(x) mutate(x, frac.abundance=(abundance/sum(as.numeric(abundance)))))
  dfs.w.p <- lapply(dfs.w.p, add.many.ps, gamma.shape = gamma.shape, gamma.rate = gamma.rate, pool.size=pool.size, n.iterations = n.iterations)
  keqs1 <- bind_rows(mutate(dfs.w.p$cele, species='cele', epsilon=median(as.matrix(cele.eps[,-1]))), 
                     mutate(dfs.w.p$human, species='human', epsilon=median(as.matrix(human.eps[,-1]), na.rm=T)), 
                     mutate(dfs.w.p$zf, species='zfish', epsilon=median(as.matrix(zf.eps[,-1])))) 
  keqs1 <- keqs1 %>%
    mutate(single.keq=p*((epsilon*(frac.abundance/pt))-(t.usage/pt))) 
  return(keqs1)
  }


# taking the pairs
add.pairs.info <- function(hex.df){
  hex.df <- hex.df %>%
    mutate(rev.temp=sapply(template, rev.comp),
           rev.prim=sapply(primer, rev.comp)) %>%
    mutate(pair=ifelse(paste0(template,'.',primer)>paste0(rev.temp,'.',rev.prim), paste0(paste0(template,'.',primer), '.',paste0(rev.temp,'.',rev.prim)), paste0(paste0(rev.temp,'.',rev.prim), '.', paste0(template,'.',primer)))) 
  return(hex.df)
}


# %>%
#   dplyr::select(template, single.keq, species) %>%
#   filter(template=='TATCTC') 
#   ggplot(., aes(single.keq, fill=species)) + geom_histogram()
  
# while (i<n.iterations) {
#   i <- i+1
#   dfs.w.p <- lapply(list(cele=cele.df, human=human.df, zf=zf.df), function(x) mutate(x, p=(sample(rgamma(pool.size,shape=gamma.shape,rate=gamma.rate)/pool.size, n()))))
#   keqs <- bind_rows(mutate(dfs.w.p$cele, species='cele'), mutate(dfs.w.p$human, species='human'), mutate(dfs.w.p$zf, species='zfish')) %>%
#     # mutate(p=(sample(rgamma(pool.size,shape=gamma.shape,rate=gamma.rate)/pool.size, n()))) %>%
#     mutate(single.keq=p*((epsilon*(abundance/pt))-(t.usage/pt))) %>%
#     dplyr::select(template, single.keq)
#   keqs1 <- inner_join(keqs1, keqs, by='template')
# }
# head(keqs1)
# group_by(template,primer) %>%
#   summarise(mean.keq = mean(single.keq), sd.keq = sd(single.keq))



# cele.df.wp <- add.many.ps(cele.df, gamma.shape = gamma.shape, gamma.rate = gamma.rate, pool.size, n.iterations = 100)
# 
# cele.eps
# 
# keqs %>%
#   inner_join(., t.groups, by='template') %>%
#   mutate(group=as.numeric(group),
#          mean.keq=-mean.keq) %>%
#   inner_join(.,melt(keq.zf, variable.name='group', value.name = 'keq.system'), by='group') %>%
#   group_by(template) %>%
#   summarise(mean.system=mean(keq.system), mean.keq=first(mean.keq), sd.system=sd(keq.system), sd.keq=first(sd.keq)) %>%
#   melt(id.vars='template') %>%
#   mutate(type=gsub(variable, pattern = '\\..+', replacement = ''),
#          method=gsub(variable, pattern = '.+.\\.', replacement = '')) %>%
#   dplyr::select(-variable) %>%
#   dcast(template+method ~ type) %>%
#   filter(template %in% sample(t.groups$template, 50)) %>%
#   ggplot(., aes(template,mean, group=method, color=method)) + geom_point() +
#   geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd))
# 
# compare.zf.cele %>%
#   filter(species=='zf') %>%
#   # filter(group %in% sample(compare.zf.cele$group, 1000)) %>%
#   rename(template=group)
#   inner_join(., keqs, by='template') %>%
#   group_by(template) %>%
#   summarise(mean.system=mean(dg), keq=first(mean.keq), sd.system=sd(dg), sd.keq=first(sd.keq))
# 
# 

# epsilon=0.05
# # pt.data <- load.pt.data(ptCounts.file = "~/mnt/edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-CeleTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv", diag.pairs = F)
# # kmer.ab.data <- load.kmer.abundance('~/mnt/edann/hexamers/genomes_kmers/WBcel235.kmerAbundance.csv')
# # deltaG.data <- load.modelled.deltaG("~/mnt/edann/hexamers/rand_hex_deltaG_ions.txt.gz")
# # 
# # matches <- filter(pt.data$matches, substr(ptPair,1,6)==substr(ptPair, 8,100))
# # hex.df <- join.pt.data(matches, pt.data$t.usage, kmer.ab.data, deltaG.data)
