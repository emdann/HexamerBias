#############################################
###### HEXAMER BINDING MODEL FUNCTIONS ######
#############################################

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(tibble)
library(ggrepel)
library(gtools)
library(RColorBrewer)
# library(fitdistrplus)

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

make_pair_df <- function(dgMat){
  pairDf <- dgMat %>% 
    melt(value.name = 'dG') %>% 
    mutate(ptPair=paste0(Var1,'.',Var2)) %>% 
    dplyr::select(ptPair, dG)
  return(pairDf)
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


compute.primer.usage <- function(pt.all.df){
  primer.usage <- pt.all.df %>%
    group_by(primer) %>%
    summarise(p.usage= sum(pt)) %>%
    rename(hex=primer)
  return(primer.usage)
  }

# tab <- data.frame(pool = factor())
# for (n in 1:100) {
#   pool <- simulate.primer.pool(10000)
#   count <- as.data.frame(table(pool))
#   tab <- full_join(tab, count, by='pool')
# }

#### Chi-SQUARE MINIMIZATION BASED ON CALORIMETRY DATA ####

epsilon.minimize.chisq <- function(pt.df, max, min=0, prob=batch.prob.uniform(), plot=T){
  min.epsilon <- pt.df %>%
    mutate(ep=t.usage/abundance) %>%
    top_n(n = 1,ep) %>%
    sample_n(1) %>%
    .$ep
  if(min < min.epsilon){min <- min.epsilon}
  best.eps.ix <- 0
  while(best.eps.ix <= 2){
    chis <- c()
    for (eps in seq(min, max, length.out=2000) ) {
      chi.sq <- rownames_to_column(data.frame(prob), var = 'primer') %>%
        rename(p=prob) %>%
        inner_join(.,pt.df, by='primer' ) %>%
        filter(pt!=0) %>%
        mutate(chi=(-dG/0.59)
               - log(4*p) 
               - log(abundance-(t.usage/eps)) 
               + log(pt/eps)
        ) %>%
        summarise(chi=sum(chi^2))
      chis <- c(chis, chi.sq)  
    }
    best.eps.ix <- which.min(unlist(chis))
    if(best.eps.ix<=2){
      max <- seq(min,max, length.out = 100)[10]  # Take 10% of scale
    }
  }
  if(plot){plot(seq(min, max, length.out=2000), unlist(chis), pch='.', xlab='epsilon', ylab='chi^2')}
  best.eps <- seq(min, max, length.out=2000)[which.min(unlist(chis))] 
  return(best.eps)
}

compute.keqs <- function(pt.df, eps, filter.pt=200){
  keqs <- pt.df %>%
    filter(pt>filter.pt) %>%
    mutate(keq=((4^6)/4)*(pt/(eps*abundance-t.usage))) 
  return(keqs)
}

predict.coverage <- function(keqs.df, eps, prob=batch.prob.uniform()){
  pred.cov <- 
    # keqs.df %>%
    rownames_to_column(data.frame(prob), var = 'primer') %>%
    rename(p=prob) %>%
    inner_join(.,keqs.df, by='primer' )%>%
    # mutate(p=prob) %>%
    mutate(phi=p*keq,
           nuc=sapply(template, prevalent_nucleotide),
           epsilon=eps) %>%
    group_by(template) %>%
    summarise(abundance=first(abundance), epsilon=first(epsilon), t.usage=first(t.usage), sum.phi=sum(phi), nuc=first(nuc)) %>%
    mutate(pred.cov=epsilon*abundance*(sum.phi/1+sum.phi)) 
  return(pred.cov)
}

plot.prediction <- function(pred.cov.df){
  # pcc.pred <- cor(pred.cov.df$t.usage, pred.cov.df$pred.cov)
  pl <- pred.cov.df %>%
    mutate(nuc=sapply(template, prevalent_nucleotide)) %>%
    ggplot(., aes(log(t.usage/sum(t.usage)), log(pred.cov/sum(pred.cov)), color=nuc)) + 
    geom_point(alpha=0.4) +
    geom_abline(slope=1, intercept=0, color='red') +
    theme_bw() +
    xlab("log(observed cov)") + ylab("log(predicted cov)") +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text = element_text(size=20),
          axis.title = element_text(size=30),
          title = element_text(size=30)) +
    geom_text(aes(x = Inf, y = -Inf, label=paste("R.sq. =",round(cor,2))), 
              hjust   = +1, 
              vjust   = -1,
              color='black') 
  return(pl)
}

#### TESTS ####

compute.keqs.noeps <- function(pt.df, filter.pt=200){
  keqs <- pt.df %>%
    filter(pt>filter.pt) %>%
    mutate(keq=((4^6)/4)*(pt/(abundance-t.usage))) 
  return(keqs)
}

predict.coverage.noeps <- function(keqs.df, prob=4/(4^6)){
  pred.cov <- keqs.df %>%
    mutate(p=prob) %>%
    mutate(phi=p*keq,
           nuc=sapply(template, prevalent_nucleotide)) %>%
    group_by(template) %>%
    summarise(abundance=first(abundance), t.usage=first(t.usage), sum.phi=sum(phi), nuc=first(nuc)) %>%
    mutate(pred.cov=abundance*(sum.phi/1+sum.phi)) 
  return(pred.cov)
}

#####################################################################
###### PRIMER POOL AND HEXAMER SEQUENCE MANIPULATION FUNCTIONS ######
#####################################################################

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


#### A BIG BUNCH OF THINGS THAT DIDN'T WORK ####

## SCALING FACTOR ESTIMATION ASSUMING DISTRIBUTION OF PRIMER CONCENTRATION ##
group.coeffs.wPrimer.all <- function(hex.df, groupsOI, gamma.shape=24.05, gamma.rate=0.98, pool.size=sum(table(pool1)), imposeP=F){
  groups.df <- hex.df %>% 
    inner_join(., t.groups, by='template') %>%
    filter(group %in% groupsOI) %>%
    mutate(dG=dG/0.59) 
  if (imposeP) {
    groups.df <- mutate(groups.df, p=4/(4^6))
  } else {
    p.conc <- data.frame(primer=t.groups$template,   p=(sample(rgamma(pool.size, 
                                                                      shape=gamma.shape, 
                                                                      rate=gamma.rate)/pool.size, 4096)))
    groups.df <- inner_join(groups.df,p.conc, by='primer')
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

epsilon.iterative <- function(pt.df, estimation.it=10, tot.it=20, sample.size=20, imposeP=F){
  i<-1
  eps.df <- data_frame(n=1:estimation.it, estimate.epsilon(pt.df, primer.pool, sample(seq(1,40),sample.size), iterations = estimation.it, imposeP=imposeP))
  colnames(eps.df)[ncol(eps.df)] <- paste0('it.',i)
  while(i<tot.it){
    i<-i+1
    eps.df <- full_join(eps.df, data_frame(n=1:estimation.it, estimate.epsilon(pt.df, primer.pool, sample(seq(1,40),sample.size), iterations = estimation.it, imposeP=imposeP)), by='n')
    colnames(eps.df)[ncol(eps.df)] <- paste0('it.',i)
  }
  return(eps.df)
}

estimate.eps.one.group <- function(hex.df, groupOI,imposeP=T){
  groups.df <- hex.df %>% 
    inner_join(., t.groups, by='template') %>%
    filter(group==groupOI) %>%
    mutate(dG=dG/0.59
           # abundance=abundance/sum(abundance)
    )  
  if (imposeP) {
    groups.df <- mutate(groups.df, p=4/(4^6))
  } else {
    groups.df <- mutate(groups.df, 
                        p=(sample(rgamma(pool.size, 
                                         shape=gamma.shape, 
                                         rate=gamma.rate)/pool.size, n())))
  }
  groups.df <- groups.df %>%
    mutate(A=sum((p^2)*(abundance^2/pt^2)),
           C=sum(((abundance*t.usage)/(pt^2))*(p^2))) %>%
    # group_by(group) %>%
    mutate(B=sum((p*abundance)/pt),
           D=sum((t.usage/pt)*p),
           N=n()
           # C.sing=sum(((abundance*t.usage)/(pt^2))*(p^2))
    ) %>%
    mutate(eps=(-D*B+N*C)/(-(B^2)+N*A)) %>%
    dplyr::select(group,eps) %>%
    summarise_all(funs(first)) 
  return(groups.df)  
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

compute.keqs.2 <- function(pt.dfs = list(cele=cele.df, human=human.df, zf=zf.df), cele.eps, human.eps, zf.eps, gamma.shape, gamma.rate, pool.size, n.iterations=10, take.pairs=F){
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

compute.keqs.fixedP <- function(pt.df, mean.eps, 
                                prob=data.frame(p=sapply(as.character(t.groups$template), primer.prob)) %>% tibble::rownames_to_column(var='primer') , 
                                # n.iterations=10, 
                                take.pairs=F){
  tot.ab <- pt.df %>% group_by(template) %>% summarise(ab=first(abundance)) %>% mutate(tot=sum(as.numeric(ab))) %>% sample_n(1) %>% .$tot 
  if (take.pairs) {
    df.w.p <- add.pairs.info(pt.df)
  } else {
    df.w.p <- pt.df
  }
  # dfs.w.p <- lapply(dfs.w.p, function(x) mutate(x, frac.abundance=(abundance/sum(as.numeric(abundance)))))
  # dfs.w.p <- lapply(dfs.w.p, add.many.ps, gamma.shape = gamma.shape, gamma.rate = gamma.rate, pool.size=pool.size, n.iterations = n.iterations)
  keqs1 <- mutate(df.w.p, 
                  epsilon=mean.eps,
                  frac.abundance=(abundance/tot.ab)
  ) %>%
    inner_join(., prob, by='primer' ) %>%
    mutate(p=4*p) %>%
    mutate(single.keq=p*((epsilon*(frac.abundance/pt))-(t.usage/pt))) 
  # mutate(single.keq=pt/(p*(epsilon*(frac.abundance)-t.usage)))
  return(keqs1)
}

compute.keqs.fixedP.totabundance <- function(pt.df, mean.eps, 
                                             prob=data.frame(p=sapply(as.character(t.groups$template), primer.prob)) %>% tibble::rownames_to_column(var='primer') , 
                                             # n.iterations=10, 
                                             take.pairs=F){
  # tot.ab <- pt.df %>% group_by(template) %>% summarise(ab=first(abundance)) %>% mutate(tot=sum(as.numeric(ab))) %>% sample_n(1) %>% .$tot 
  if (take.pairs) {
    df.w.p <- add.pairs.info(pt.df)
  } else {
    df.w.p <- pt.df
  }
  # dfs.w.p <- lapply(dfs.w.p, function(x) mutate(x, frac.abundance=(abundance/sum(as.numeric(abundance)))))
  # dfs.w.p <- lapply(dfs.w.p, add.many.ps, gamma.shape = gamma.shape, gamma.rate = gamma.rate, pool.size=pool.size, n.iterations = n.iterations)
  keqs1 <- mutate(df.w.p, 
                  epsilon=mean.eps
  ) %>%
    inner_join(., prob, by='primer' ) %>%
    mutate(p=4*p) %>%
    mutate(single.keq=p*((epsilon*(abundance/pt))-(t.usage/pt))) 
  # mutate(single.keq=pt/(p*(epsilon*(frac.abundance)-t.usage)))
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

