### RUN COVERAGE PREDICTION ### 
library(purrr)
library(tidyr)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("inputPattern", type="character",
                    help = "Pattern of RDS files to be processed")
parser$add_argument("-o", "--outputPrefix", type="character",
                    help="prefix for output file")
args <- parser$parse_args()

setwd("/hpc/hub_oudenaarden/edann")
input.dir <- "./VAN2591/"

input.pattern <- args$inputPattern
output.prefix <- args$outputPrefix

keqs.file <- "./mean_keqs_noBS_all.RData"
eps.model.file <- "./model_epsilon.RData"
gencov.file <- "./VAN2591/genomecov_all.txt"

rdata.pattern <- paste0(input.pattern, '.+.RDS')
input.files <- list.files(input.dir, pattern = rdata.pattern, full.names = T)
load(eps.model.file)
load(keqs.file)

pt.R1.all.df <- readRDS(input.files[grep(input.files, pattern = 'handMixOld')])
pt.R2.all.df <- readRDS(input.files[grep(input.files, pattern = 'handMixNew')])
pt.T.all.df <- readRDS(input.files[grep(input.files, pattern = 'moreT')])
pt.G.all.df <- readRDS(input.files[grep(input.files, pattern = 'moreG')])

gencov.noBS <- read.table(gencov.file, sep='\t', header=F, col.names = c('smp', 'genomecov'))
gencov.noBS <- gencov.noBS %>%
  mutate(genomecov=1-genomecov) %>%
  filter(grepl(smp, pattern = 'noBS')) 

eps.pt.handMixOld <- exp(predict(model.epsilon, data.frame(genomecov=filter(gencov.noBS, grepl(smp, pattern = paste0(input.pattern, '.+handMixOld')))$genomecov)))
eps.pt.handMixNew <- exp(predict(model.epsilon, data.frame(genomecov=filter(gencov.noBS, grepl(smp, pattern = paste0(input.pattern, '.+handMixNew')))$genomecov)))
eps.pt.moreG <- exp(predict(model.epsilon, data.frame(genomecov=filter(gencov.noBS, grepl(smp, pattern = paste0(input.pattern, '.+moreG')))$genomecov)))
eps.pt.moreT <- exp(predict(model.epsilon, data.frame(genomecov=filter(gencov.noBS, grepl(smp, pattern = paste0(input.pattern, '.+moreT')))$genomecov)))

pt.R1.keqs <- inner_join(select(avg.keqs.all, primer, template, keq), pt.R1.all.df, by=c("primer", 'template')) 
pt.R2.keqs <- inner_join(select(avg.keqs.all, primer, template, keq), pt.R2.all.df, by=c("primer", 'template')) 
pt.G.keqs <- inner_join(select(avg.keqs.all, primer, template, keq), pt.G.all.df, by=c("primer", 'template')) 
pt.T.keqs <- inner_join(select(avg.keqs.all, primer, template, keq), pt.T.all.df, by=c("primer", 'template'))

prob.t <- batch.prob.uniform(nuc.probs = c(pA=0.25, pT=0.45, pC=0.25, pG=0.05))
prob.g <- batch.prob.uniform(nuc.probs = c(pA=0.25, pG=0.45, pC=0.25, pT=0.05))

pred.cov.pt.G <- predict.coverage(pt.G.keqs, eps = eps.pt.moreG, prob = prob.g)
pred.cov.pt.T <- predict.coverage(pt.T.keqs, eps = eps.pt.moreT, prob = prob.t)
pred.cov.pt.R1 <- predict.coverage(pt.R1.keqs, eps = eps.pt.handMixOld)
pred.cov.pt.R2 <- predict.coverage(pt.R2.keqs, eps = eps.pt.handMixNew)

pred.coverage <- list(handMixOld=pred.cov.pt.R1,
                      handMixNew=pred.cov.pt.R2,
                      moreG=pred.cov.pt.G,
                      moreT=pred.cov.pt.T
                      )

### R-sq w different primer batch 
g.probs <- seq(0,0.5,0.05)
prob.vecs <- lapply(g.probs, function(x) batch.prob.uniform(nuc.probs = c(pA=0.25, pT=1-0.5-x, pG=x, pC=0.25)))

variable.batch.Rsq <- function(pt.keqs, eps){
  cor.pt <- map(prob.vecs, predict.coverage, keqs.df=pt.keqs, eps = eps) %>%
    map_dbl(function(x) cor(x$pred.cov, x$t.usage))
  return(cor.pt)
  }

cor.pt.R1 <- variable.batch.Rsq(pt.R1.keqs, eps.pt.handMixOld)
cor.pt.R2 <- variable.batch.Rsq(pt.R2.keqs, eps.pt.handMixNew)
cor.pt.G <- variable.batch.Rsq(pt.G.keqs, eps.pt.moreG)
cor.pt.T <- variable.batch.Rsq(pt.T.keqs, eps.pt.moreT)

pcc.primer.batch <- data.frame(prob.G = g.probs, prob.T = 1-g.probs, 
           moreT = cor.pt.T, 
           moreG = cor.pt.G,
           handMixOld = cor.pt.R1,
           handMixNew = cor.pt.R2) %>%
  gather(key='sample', value='PCC', 3:5)
  
output <- list(pred.coverage=pred.coverage, 
               pcc.primer.batch=pcc.primer.batch)
save(output,
  file=paste0(input.dir, output.prefix, '.primerbatch.predcoverage.RData'))









