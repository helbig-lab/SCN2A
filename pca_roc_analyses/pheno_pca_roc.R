library(tidyverse)
library(Hmisc)
library(logisticPCA)
library(ROCR)
library(corrplot)
library(RColorBrewer)

start <- Sys.time()
message(" \n Starting PCA pheno anaylsis... \n ")

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
# Load pca
################## 

# Output from 'run_pheno_ks.R'
load("pca_roc_analyses/phenogroup_logpca_model_k3.RData")

################## 
# Analysis
################## 

pats = as.vector(scn2a_prop$famID %>% unique())
hpo = as.vector(scn2a_prop$HPO %>% unique())

# PC values
pca_x <- as.data.frame(logpca_model$PCs); pca_x$famID = pats
colnames(pca_x) = c("PC1", "PC2", "PC3", "famID")
pca_x <- pca_x %>% select(famID, 1:(ncol(pca_x)-1)) %>% left_join(scn2a_prop %>% select(-HPO) %>% unique())

# Variance explained
logpca_model$prop_deviance_expl # total variance
vars_transformed <- apply(pca_x %>% select("PC1", "PC2", "PC3"), 2, var)
vars_transformed/sum(vars_transformed)*logpca_model$prop_deviance_expl

# Factor loadings
loadings <- as.data.frame(logpca_model$U); loadings$HPO = hpo
colnames(loadings) = c("PC1", "PC2", "PC3", "HPO")
loadings <- loadings %>% left_join(hpo_all) %>% select(HPO, def, 1:(ncol(loadings)-1))

###

start <- Sys.time()
message(" \n Starting ROC pheno anaylsis... \n ")

####################################
# 2. ROC
####################################

#########
# Data input
#########

# Requires PC values from PCA

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
# ROC phenotypic group binary classification
#########

roc_x_1 <- roc_x %>%
  mutate(broad_phx = case_when(broad_phx == "benign" ~ "BFNIS",
                               broad_phx == "encephalopathy" ~ "DEE",
                               broad_phx == "ASD" ~ "Autism"))

g <- c("DEE", "Autism", "BFNIS")
summary_roc <- expand.grid(PHENO = g, PC = 1:3)
summary_roc$SENSITIVITY = NA
summary_roc$SPECIFICITY = NA
summary_roc$OPT_PC_CUTOFF = NA
summary_roc$AUC = NA

for (pcx in unique(summary_roc$PC)) {
  
  for (i in 1:length(g)) {
    
    g1 <- g[i]
    
    # Predict groups as 1 or 0
    df <- roc_x_1 %>% select(famID, paste0("PC", pcx), broad_phx)
    
    m1 = mean(df %>% filter(broad_phx == g1) %>% select(paste0("PC", pcx)) %>% pull())
    m2 = mean(df %>% filter(broad_phx %nin% g1) %>% select(paste0("PC", pcx)) %>% pull())
    
    if (m1 < m2) {
      df <- df %>% mutate(pred = case_when(broad_phx == g1 ~ 0, TRUE ~ 1))
    } else {
      df <- df %>% mutate(pred = case_when(broad_phx == g1 ~ 1, TRUE ~ 0))
    }
    
    pred_std <- prediction(df %>% select(paste0("PC", pcx)), df$pred)
    perf_std <- performance(pred_std,"tpr","fpr")
    
    cutoffs_std <- data.frame(cutoff=perf_std@alpha.values[[1]], fpr=perf_std@x.values[[1]], tpr=perf_std@y.values[[1]])
    opt_std = opt.cut(perf_std, pred_std)
    sen=opt_std[1]
    spec=opt_std[2]
    cut = opt_std[3]
    
    # Get the area under the curve (AUC)
    auc.perf_std = performance(pred_std, measure = "auc")
    auc = auc.perf_std@y.values[[1]]

    summary_roc[summary_roc$PHENO == g1 & summary_roc$PC == pcx,]$SENSITIVITY = sen
    summary_roc[summary_roc$PHENO == g1 & summary_roc$PC == pcx,]$SPECIFICITY = spec
    summary_roc[summary_roc$PHENO == g1 & summary_roc$PC == pcx,]$AUC = auc
    summary_roc[summary_roc$PHENO == g1 & summary_roc$PC == pcx,]$OPT_PC_CUTOFF = cut
    
  }
  
}

# Output for ROC performance for phenogroups
write.csv(summary_roc, paste0(input.yaml$output_dir,"pca_roc_analyses/phenogroup_roc_summary.csv"))

#########
# ROC phenotypic group specific comparisons
#########

roc_x_2 <- roc_x %>% 
  filter(broad_phx %in% c("benign", "encephalopathy", "ASD")) %>% 
  mutate(broad_phx = case_when(broad_phx == "benign" ~ "BFNIS",
                               broad_phx == "encephalopathy" ~ "DEE",
                               broad_phx == "ASD" ~ "Autism"))

for (pcx in 1:3) {
  
  corr_table <- as.data.frame(matrix(NA,nrow=length(roc_x_2$broad_phx %>% unique()),ncol=length(roc_x_2$broad_phx %>% unique())))
  row.names(corr_table) = roc_x_2$broad_phx %>% unique()
  colnames(corr_table) = roc_x_2$broad_phx %>% unique()
  
  for (i in 1:(NROW(corr_table)-1)) {
    for (j in (i+1):NCOL(corr_table)) {
      
      g1 = row.names(corr_table)[i]
      g2 = colnames(corr_table)[j]
      
      if (g1 == g2) {
        auc = 0
      } else {
        
        # Predict groups as 1 or 0
        df <- roc_x_2 %>% 
          select(famID, paste0("PC", pcx), broad_phx) %>% 
          filter(broad_phx %in% c(g1,g2)) 
        
        m1 = mean(df %>% filter(broad_phx == g1) %>% select(paste0("PC", pcx)) %>% pull())
        m2 = mean(df %>% filter(broad_phx == g2) %>% select(paste0("PC", pcx)) %>% pull())
        
        if (m1 < m2) {
          df <- df %>% mutate(pred = case_when(broad_phx == g1 ~ 0, broad_phx == g2 ~ 1))
        } else {
          df <- df %>% mutate(pred = case_when(broad_phx == g1 ~ 1, broad_phx == g2 ~ 0))
        }
        
        pred_std <- prediction(df %>% select(paste0("PC", pcx)), df$pred)
        perf_std <- performance(pred_std,"tpr","fpr")
        
        cutoffs_std <- data.frame(cutoff=perf_std@alpha.values[[1]], fpr=perf_std@x.values[[1]], tpr=perf_std@y.values[[1]])
        opt_std = opt.cut(perf_std, pred_std)
        sen=opt_std[1]
        spec=opt_std[2]
        
        # Get the area under the curve (AUC)
        auc.perf_std = performance(pred_std, measure = "auc")
        auc = auc.perf_std@y.values[[1]]
        
      }
      
      corr_table[i,j] = auc
      
    }
  }
  
  corr_table[is.na(corr_table)] = 0
  write.csv(corr_table, paste0(input.yaml$output_dir,"pca_roc_analyses/pheno_roc_performance_PC", pcx, "_matrix.csv"), row.names = T)
  
}

message("\n  ...pheno PCA-ROC analysis complete \n ")
stop = Sys.time()
stop - start
