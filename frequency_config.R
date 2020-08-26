library(tidyverse)
library(memoise)

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
  message('\n  Please mention the HPO Tree File (Should be able to get it from [here]("https://github.research.chop.edu/KAUFMANMC/SCN2A/tree/master/raw_files")  ) - Cant Proceed without that \n')
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
  variant_file <- read_csv(input.yaml$variant_file)
}else{
  message('\n  Please mention the Variants File in the specified format (Should be able to get it from [here]("https://github.com/helbig-lab/SCN2A/tree/master/raw_files") ) - Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$freq_dir) == T){
  message('\n Please determine frequency analysis directory in input config file - Cant Proceed without that \n')
  break;
