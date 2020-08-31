library(tidyverse)
library(Hmisc)
library(logisticPCA)
library(reshape2)
library(ROCR)
library(corrplot)
library(RColorBrewer)

start <- Sys.time()
message(" \n Starting PCA variant anaylsis... \n ")

setwd(input.yaml$pca_dir)

####################################
# 1. PCA
####################################

################## 
# Data input
################## 

scn2a <- read.csv(input.yaml$variant_file, stringsAsFactors = F) %>% 
  filter(variant_type_1 != "exclude") %>%
  filter(grepl("HP:[0-9]", HPO)) %>% 
  unique()

func_variants <- read.csv(input.yaml$functional_variant, stringsAsFactors = F) %>% 
  rename(effect = Overall.Effect) %>% 
  select(variant, effect)

scn2a_pos <- read.csv(input.yamlpos_prop, stringsAsFactors = FALSE)
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

scn2a <- scn2a %>% filter(!is.na(variant_type_2))


################## 
# Load pca
################## 

# Output from 'run_variant_ks.R'
load("pca_roc_analyses/functional_variant_logpca_model_k3.RData")

################## 
# Analysis
################## 

pats = as.vector(scn2a$famID %>% unique())
hpo = as.vector(scn2a_prop$HPO %>% unique())

# PC values
pca_x <- as.data.frame(logpca_model$PCs); pca_x$famID = pats
colnames(pca_x) = c("PC1", "PC2", "PC3", "famID")
pca_x <- pca_x %>% select(famID, 1:(ncol(pca_x)-1)) %>% left_join(scn2a_prop %>% select(-HPO) %>% unique())

# Variance explained
logpca_model$prop_deviance_expl
vars_transformed <- apply(pca_x %>% select("PC1", "PC2", "PC3"), 2, var)
vars_transformed/sum(vars_transformed)*logpca_model$prop_deviance_expl

# Factor loadings
loadings <- as.data.frame(logpca_model$U); loadings$HPO = hpo
colnames(loadings) = c("PC1", "PC2", "PC3", "HPO")
loadings <- loadings %>% left_join(hpo_all) %>% select(HPO, def, 1:(ncol(loadings)-1))

###

start <- Sys.time()
message(" \n Starting ROC variant anaylsis... \n ")

####################################
# 2. ROC
####################################

#########
# Data input
#########

# Requires PC values

roc_x <- pca_x
# rm(list=setdiff(ls(), "roc_x"))

# Get "optimal" cut point
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

#########
# ROC and PPV analysis
#########

roc <- list()
auc_all <- c()
for (pcx in 1:3) {
  
  g1 = "Missense LOF + PTV"
  g2 = "Missense GOF"
  
  # Decide what group to use 0 and 1
  df <- roc_x %>% 
    select(famID, paste0("PC", pcx), functional_1) %>% 
    filter(functional_1 %in% c(g1, g2)) 
  
  if (pcx == 2) {
    df_PC2 = df
  }
  
  m1 = median(df %>% filter(functional_1 == g1) %>% select(paste0("PC", pcx)) %>% pull())
  m2 = median(df %>% filter(functional_1 == g2) %>% select(paste0("PC", pcx)) %>% pull())
  
  if (m1 < m2) {
    df <- df %>% mutate(pred = case_when(functional_1 == g1 ~ 0, functional_1 == g2 ~ 1))
  } else {
    df <- df %>% mutate(pred = case_when(functional_1 == g1 ~ 1, functional_1 == g2 ~ 0))
  }
  
  pred_std <- prediction(df %>% select(paste0("PC", pcx)), df$pred)
  perf_std <- performance(pred_std,"tpr","fpr")
  cutoffs_std <- data.frame(cutoff=perf_std@alpha.values[[1]], fpr=perf_std@x.values[[1]], tpr=perf_std@y.values[[1]])
  
  # Calculate PPV for PC2
  if (pcx == 2) {
    
    ppv_std <- cutoffs_std
    ppv_std$ppv_gof = NA
    ppv_std$ppv_lof = NA
    for (c in 1:(nrow(ppv_std))) {
      
      cut <- ppv_std$cutoff[c]
      
      cut_l <- df %>% filter(PC2 <= cut)
      cut_r <- df %>% filter(PC2 > cut)
      
      if (nrow(cut_r) > 0 & nrow(cut_l) > 0) {
        ppv_std$ppv_gof[c] <- sum(cut_r == 1)/(sum(cut_r == 1) + sum(cut_r == 0))
        ppv_std$ppv_lof[c] <- sum(cut_l == 0)/(sum(cut_l == 1) + sum(cut_l == 0))
      } else if (nrow(cut_l) == 0) {
        ppv_std$ppv_gof[c] <- sum(cut_r == 1)/(sum(cut_r == 1) + sum(cut_r == 0))
        ppv_std$ppv_lof[c] <- NA
      } else if (nrow(cut_r) == 0) {
        ppv_std$ppv_lof[c] <- sum(cut_l == 0)/(sum(cut_l == 1) + sum(cut_l == 0))
        ppv_std$ppv_gof[c] <- NA
      }
      
    }
    ppv_std <- ppv_std %>% na.exclude()
    
    ppv_df <- ppv_std %>%
      select(cutoff, ppv_lof) %>%
      rename(ppv = ppv_lof) %>%
      mutate(pred = "Missense LOF + PTV") %>%
      rbind(ppv_std %>%
              select(cutoff, ppv_gof) %>%
              rename(ppv = ppv_gof) %>%
              mutate(pred = "Missense GOF"))
  }
  
  cutoffs_std <- cutoffs_std %>% select(-cutoff)
  roc[[pcx]] <- cutoffs_std
  
  opt_std = opt.cut(perf_std, pred_std)
  sen=opt_std[1]
  spec=opt_std[2]
  
  # Get the area under the curve (AUC)
  auc.perf_std = performance(pred_std, measure = "auc")
  auc = auc.perf_std@y.values[[1]]
  auc_all <- append(auc_all, auc)
  
}

# ROC plot (PC 1-3)

df.long <- melt(roc,id.vars ="fpr") %>% rename(PC = L1) %>% mutate(PC = paste0("PC",PC))

out_roc <- df.long %>% 
  mutate(tpr = value) %>% 
  select(PC, tpr, fpr) %>% 
  mutate(AUC = case_when(PC == "PC1" ~ auc_all[1],
                         PC == "PC2" ~ auc_all[2],
                         PC == "PC3" ~ auc_all[3]))

write.csv(out_roc, paste0(input.yaml$output_dir,"pca_roc_analyses/variant_ROC.csv"))

###
# PPV data for PC2
###

variants <- roc_x %>% select(famID, PC2, variant, variant_type_2, functional_1) %>% 
  left_join(scn2a %>% select(famID, broad_phx) %>% unique())

variants$cutoff = NA
for (i in 1:NROW(variants)) {
  variants$cutoff[i] <- ppv_df[which(abs(ppv_df$cutoff-variants$PC2[i]) == min(abs(ppv_df$cutoff-variants$PC2[i]))),]$cutoff %>% unique()
}
ppv_variants <- variants %>% left_join(ppv_df)

## Table
ppv_variants_1 <- variants %>% 
  left_join(ppv_std) %>% 
  select(famID, variant, functional_1, PC2, ppv_gof, ppv_lof)

count <- ppv_variants_1 %>% 
  dplyr::group_by(variant) %>% 
  count()

ppv_variants_1 <- ppv_variants_1 %>% 
  left_join(count) %>% 
  rename(variant_n = n) %>% 
  select(1,2,7,3,4,5,6)

ppv_variants_1[ppv_variants_1$variant == "E1211K",]$functional_1 = "Mixed"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense GOF",]$functional_1 = "GoF"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense LOF + PTV",]$functional_1 = "LoF"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense other",]$functional_1 = "none"

write.csv(ppv_variants_1, paste0(input.yaml$output_dir,"pca_roc_analyses/variant_pred_PC2_PPV.csv", row.names = F))

message("\n  ...variant PCA-ROC analysis complete \n ")
stop = Sys.time()
stop - start

