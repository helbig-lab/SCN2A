library(tidyverse)
library(memoise)

if(is.null(input.yaml$functional_variant) == T){
  message('\n Please determine functional variant file in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$pos_prop) == T){
  message('\n Please determine positive propagation file in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$pos_prop) == T){
  message('\n Please determine negative propagation file in input config file - Cant Proceed without that \n')
  break;
}

system(paste0("mkdir -p ",input.yaml$output_dir,"pca_roc_analyses"))
system(paste0("cp scripts/run_pheno_k3.R ",input.yaml$output_dir,"pca_roc_analyses/."))
system(paste0("cp scripts/run_variant_k3.R ",input.yaml$output_dir,"pca_roc_analyses/."))
system(paste0("cp scripts/pheno_pca_roc.R ",input.yaml$output_dir,"pca_roc_analyses/."))
system(paste0("cp scripts/variant_pca_roc.R ",input.yaml$output_dir,"pca_roc_analyses/."))
setwd(paste0(input.yaml$output_dir,"pca_roc_analyses"))

  source(paste0("run_pheno_k3.R"))
  source(paste0("run_variant_k3.R"))
  source(paste0("pheno_pca_roc.R"))
  source(paste0("variant_pca_roc.R")) 


