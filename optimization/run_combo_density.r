library(parallel)
source("/hpc/hub_oudenaarden/edann/bin/coverage_bias/deltaGprediction/binding_model_functions.r")

compute.enrichment.score <- function(dens.df, fc.df){
  # dens.df <- pred.cov.df %>%
  #   mutate(binding.dens = pred.cov/abundance) 
  enrichment.score <- dens.df %>%
    inner_join(., fc, by="template") %>%
    mutate(dens.fc=binding.dens*fc) %>%
    summarise(score=sum(dens.fc)) %>%
    .$score
  return(enrichment.score)
}

joining.fun <- function(...){
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  col1 = colnames(df1)[1]
  col2 = colnames(df2)[1]
  xxx = left_join(..., by = setNames(col2,col1))
  return(xxx)
}

density.combo <- function(prob.vec, keqs.df=d3r.keqs, eps=epsilon.d3r){
  b.probs <- batch.prob.uniform(hexs=all.hexamers(), nuc.probs = prob.vec)
  pred.cov.b <- predict.coverage(keqs.df, eps, prob = b.probs)
  dens.df <- pred.cov.b %>%
    mutate(binding.dens = pred.cov/abundance)  %>% 
    dplyr::select(template, binding.dens)
  colnames(dens.df)[2] <- do.call(paste,c('dens',prob.vec, sep='_'))
  return(dens.df)
  }

## Load Keqs and epsilon for mouse BS-seq
load("/hpc/hub_oudenaarden/edann/primer_conc_VAN2493/d3r_keqs.RData")
epsilon.d3r <- 612.6353

## Make table of sequence combosition combos
prob.combos <- hexamerMatrix(stepSize = 0.1)
test.combos <- prob.combos[which(prob.combos['pC']!=0 & prob.combos['pG']!=0),]
l.test.combos <- lapply(seq_len(nrow(test.combos)), function(i) test.combos[i,])

## Compute density for all combos and save in table
# test.combo.density <- mclapply(sample(l.test.combos, 5), density.combo, mc.cores = detectCores())
test.combo.density <- lapply(sample(l.test.combos, 5), density.combo)
save(test.combo.density, file="/hpc/hub_oudenaarden/edann/test_primer_combos_density.RData")

# dens.table <- Reduce( joining.fun, test.combo.density)
# save(dens.table, file="/hpc/hub_oudenaarden/edann/primer_combos_density.RData")
