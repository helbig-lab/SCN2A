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


if(is.null(input.yaml$output_dir) == T){
  message('\n Please mention the Field output_dir in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$file_path) == T){
  message('\n Please mention the Field file_path in input config file - Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$hpo_tree) == F){
  hpo_tree <- read_csv(input.yaml$hpo_tree)
}else{
  message('\n  Please mention the HPO Tree File (Should be able to get it from [here]("https://github.com/helbig-lab/SCN2A/tree/master/raw_files")  ) - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$hpo_ancestor) == F){
  hpo_ancestor <- read_csv(input.yaml$hpo_ancestor)
}else{
  message('\n  Please mention the HPO Ancestor File (Should be able to get it from [here]("https://github.com/helbig-lab/SCN2A/tree/master/raw_files") ) - Cant Proceed without that \n')
  break;
}

options(stringsAsFactors = F)

if(is.null(input.yaml$hpo_path) == F){
  hpo_path <- read_csv(input.yaml$hpo_path)
}else{
  message('\n  Please mention the HPO Path File (Should be able to get it from [here]("https://github.com/helbig-lab/SCN2A/tree/master/raw_files")  ) - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$variant_file) == F){
  variant <- read_csv(input.yaml$variant_file)
}else{
  message('\n  Please mention the Variants File in the specified format (Should be able to get it from [here]("https://github.com/helbig-lab/SCN2A/tree/master/raw_files") ) - Cant Proceed without that \n')
  break;
}

#Term propagation
if(is.null(input.yaml$pos_ic) == F ){
  pos_ic <- read_csv(input.yaml$pos_ic)
}
else{
    message("\n  Manually propagating information content... \n ")
    source("raw_files/compose_base_prop_ic.R")
    input.yaml$pos_ic <- read_csv("raw_files/pos_IC_manual.csv")
}


#Similarity analyses
if(is.null(input.yaml$sim_dir) == F ){
  message("\n  Running similarity analysis... \n ")
  source(sim_config.R)
}
else{
    message("\n  Checking for similarity analysis directory source... \n ")
    source(sim_config.R)
    if(is.null(input.yaml$sim_dir) == F){
      message("\n  Running similarity analyses... \n ")
      source(sim_config.R)
    }
    else {
      next;
    }
}

#PCA analyses
if(is.null(input.yaml$pca_dir) == F ){
  message("\n  Running PCA analysis... \n ")
  source(pca_config.R)
}
else {
  next;
}

#Frequency analyses
if(is.null(input.yaml$freq_dir) == F){
message("\n  Running frequency analyses... \n ")
source(frequency_config.R)
}
else {
  next;
}


message("\n  ...all selected analyses complete \n ")
stop = Sys.time()
stop - start


