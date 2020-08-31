library(tidyverse)
library(memoise)

start <- Sys.time()
message(" \n Starting config file... \n ")


message("\n  ...Checking for analyses to run... \n ")

#General yaml files needed for most analyses
capture <- commandArgs(trailingOnly = TRUE)

opt1 = list(make_option(c("--input"), type = "character", default = "input.yml", dest = "input"))


user_input <- function(name, argv) {
  return(name %in% names(argv))
}

argv <- parse_args(OptionParser(option_list = opt1))

if (user_input("input", argv)) {
  input = argv$input 
    if(file.exists(input) == T){
      input.yaml <- yaml::read_yaml(input)
    }else{
      message('\n Input YAML not found \n')
      break;
    }
  } else {
  message('Cannot proceed without input yaml file. Please use "--input" flag .\n')
}

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
  

if(is.null(input.yaml$pca_dir) == T){
  message('\n Please determine directory in input config file - Cant Proceed without that \n')
  break;
}
else{
  source(paste0(input.yaml$pca_dir, run_pheno_k3.R))
  source(paste0(input.yaml$pca_dir, run_variant_k3.R))
  source(paste0(input.yaml$pca_dir, pheno_pca_roc.R))
  source(paste0(input.yaml$pca_dir, variant_pca_roc.R)) 
}
