setwd("/hpc/hub_oudenaarden/")
source("./edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

tabDg <- read.delim(gzfile("./edann/hexamers/rand_hex_deltaG_ions.txt.gz"), sep=' ', header = FALSE, col.names=c('template', 'dG'))
t.groups <- read.table("./edann/MatchingTemplate_groups.tsv", col.names = c('template', 'group'))

kmer.ab.cele <- load.kmer.abundance('./edann/hexamers/genomes_kmers/WBcel235.kmerAbundance.csv')
kmer.ab.zf <- load.kmer.abundance('./edann/hexamers/genomes_kmers/danRer10.kmerAbundance.csv')
kmer.ab.human <- load.kmer.abundance('./edann/hexamers/genomes_kmers/hg38.kmerAbundance.csv')

## Loading data
print("Loading cele..")
cele.all <- load.pt.data("./edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-CeleTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv")
cele.all.df <- join.pt.data(cele.all$matches, cele.all$t.usage, kmer.ab.cele, tabDg)

print("Loading zf..")
zf.all <- load.pt.data("./edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-zfishTotal-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv")
zf.all.df <- join.pt.data(zf.all$matches, zf.all$t.usage, kmer.ab.zf, tabDg)

print("Loading human..")
human.all <- load.pt.data("./edann/VAN2423_onePreamp/cov_prediction/CG-pbat-gDNA-humanAPKS-noBS-1preAmp-handMix_lmerged_R1.ptCounts.qualFilt.csv")
human.all.df <- join.pt.data(human.all$matches, human.all$t.usage, kmer.ab.human, tabDg)

## Simulating primer pool
primer.pool <- simulate.primer.pool(pool.size = 200000)
fit.gamma <- summary(fitdist(as.numeric(table(primer.pool)), 'gamma'))
gamma.shape = fit.gamma$estimate['shape']
gamma.rate = fit.gamma$estimate['rate']
pool.size = length(primer.pool)

## Estimation of scaling factor
cele.eps <-epsilon.iterative(cele.df, estimation.it = 20)
zf.eps <-epsilon.iterative(zf.df, estimation.it = 20)
human.eps <-epsilon.iterative(human.df, estimation.it = 20)

## Filter out low counts
cele.filt.df <- filter(cele.all.df, pt > 200)
human.filt.df <- filter(human.all.df, pt > 200)
zf.filt.df <- filter(zf.all.df, pt > 200)

## Estimation of equilibrium constants for all pairs

keqs <- compute.keqs(pt.dfs=list(cele=cele.filt.df, human=human.filt.df, zf=zf.filt.df),
             cele.eps = cele.eps, human.eps = human.eps, zf.eps = zf.eps,
             gamma.shape = gamma.shape, gamma.rate = gamma.rate, pool.size = pool.size,
             n.iterations = 50, take.pairs = T)


keqs.mean <- keqs %>% 
  # filter(species=='cele') %>% 
  filter(pt>300) %>%
  group_by(pair) %>%
  summarise(dG=first(dG), mean.keq=mean(single.keq, na.rm=T), sd.keq=sd(single.keq, na.rm=T))

write.csv(x = keqs.mean, "/hpc/hub_oudenaarden/edann/Keqs_all_pairs.csv")