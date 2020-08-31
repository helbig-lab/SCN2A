# Logistic PCA for phenogroups

library(tidyverse)
library(Hmisc)
library(logisticPCA)

setwd(input.yaml$pca_dir)

################## 
# Data input
################## 

scn2a <- read.csv(input.yaml$variant_file, stringsAsFactors = F) %>% 
  filter(variant_type_1 != "exclude") %>%
  filter(grepl("HP:[0-9]", HPO)) %>% 
  unique()

scn2a_pos <- read.csv(input.yaml$pos_prop, stringsAsFactors = FALSE)
scn2a_neg <- read.csv(input.yaml$neg_prop, stringsAsFactors = FALSE)
scn2a_prop <- scn2a_pos %>% rbind(scn2a_neg)
rm(scn2a_pos, scn2a_neg)

# Hpo dictionary
pos_hpo <- read.csv(input.yaml$hpo_ancestor, stringsAsFactors = FALSE) %>% 
  select(term, def) %>% 
  rename(HPO = term) %>% 
  mutate(HPO = trimws(HPO))

################## 
# Data cleaning
################## 

scn2a_prop <- scn2a_prop %>% 
  left_join(scn2a %>% select(famID, broad_phx) %>% unique())

# Add negatives to hpo dictionary
neg_hpo <- pos_hpo %>% 
  mutate(HPO = gsub("H","N",HPO)) %>% 
  mutate(def = paste("No ",def, sep = "")) %>% 
  unique()

hpo_all <- pos_hpo %>% rbind(neg_hpo)
rm(neg_hpo)

################## 
# Set up matrix
################## 

ds <- scn2a_prop
pats = as.vector(ds$famID %>% unique())
hpo = as.vector(ds$HPO %>% unique())

mat <- data.frame(matrix(NA,nrow=length(hpo),ncol=length(pats)))
dimnames(mat)[[1]] = hpo
dimnames(mat)[[2]] = pats

for (i in 1:NROW(ds)) {
  
  p <- as.character(ds$famID[i])
  hpo <- as.character(ds$HPO[i])
  mat[rownames(mat) == hpo, colnames(mat) == p] <- 1
  
}
mat[is.na(mat)] <- 0
mat <- t(mat)

################## 
# Logistic PCA
##################

j = 3 # fit for 3 dimensions
  
logpca_cv = cv.lpca(mat, ks = j, ms = 1:10) 

# Cross validation- resampling to find optimal 'm'
logpca_model = logisticPCA(mat, k = j, m = which.min(logpca_cv)) 
  
save(logpca_model, file=paste0("phenogroup_logpca_model_k", j, ".RData"))


