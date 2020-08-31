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
