library(tidyverse)
library(Hmisc)
library(logisticPCA)

options(stringsAsFactors = F)

start <- Sys.time()
message(" \n Generating pca pheno PCA-ROC data... \n ")


################## 
# Data input
################## 

scn2a <- read.csv(input.yaml$variant_file) %>% 
  filter(variant_type_1 != "exclude") %>%
  filter(grepl("HP:[0-9]", HPO)) %>% 
  unique()

func_variants <- read.csv(input.yaml$functional_variant, stringsAsFactors = F) %>% 
  rename(effect = Overall.Effect) %>% 
  select(variant, effect)

scn2a_pos <- read.csv(input.yaml$pos_prop)
scn2a_neg <- read.csv(input.yaml$neg_prop)
scn2a_prop <- scn2a_pos %>% rbind(scn2a_neg)
rm(scn2a_pos, scn2a_neg)

# Hpo dictionary
pos_hpo <- read.csv(input.yaml$hpo_ancestor) %>% 
  select(term, def) %>% 
  rename(HPO = term) %>% 
  mutate(HPO = trimws(HPO))

################## 
# Data cleaning
################## 

variants_all <- scn2a %>% 
  select(variant, variant_type_2) %>% 
  unique() %>% 
  filter(variant_type_2 != "") %>% 
  left_join(func_variants) %>% 
  mutate(functional = case_when(variant_type_2 == "missense" & effect == "GOF" ~ "Missense GOF",
                                variant_type_2 == "missense" & effect == "LOF" ~ "Missense LOF",
                                variant_type_2 == "missense" & effect %nin% c("GOF", "LOF") ~ "Missense other",
                                variant_type_2 == "PTV" ~ "PTV")) %>% 
  mutate(functional_1 = case_when(functional == "Missense GOF" ~ "Missense GOF",
                                  functional %in% c("Missense LOF", "PTV") ~ "Missense LOF + PTV",
                                  functional == "Missense other" ~ "Missense other"))

scn2a_prop <- scn2a_prop %>% 
  left_join(scn2a %>% select(famID, variant) %>% unique()) %>% 
  left_join(variants_all %>% select(variant, variant_type_2, functional_1) %>% unique())


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

ds <- scn2a_prop %>% filter(!is.na(variant_type_2))
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

j = 3

logpca_cv = cv.lpca(mat, ks = j, ms = 1:10)
logpca_model = logisticPCA(mat, k = j, m = which.min(logpca_cv))

save(logpca_model, file=paste0("functional_variant_logpca_model_k", j, ".RData"))

message(" \n ...variant data generated \n ")
stop = Sys.time()
stop - start

