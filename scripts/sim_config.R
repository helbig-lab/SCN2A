library(tidyverse)
library(memoise)

  input.yaml$sim_cluster <- as.character(paste0(input.yaml$file_path, "sim_analyses/"))

if(input.yaml$sim_algorithm == 1){
  input.yaml$sim_dir <- as.character("sim_analyses/resnik/")
  message('\n  Similarity analysis algorithm selected - Resnik \n')
}
if(input.yaml$sim_algorithm == 2){
 input.yaml$sim_dir <- as.character("sim_analyses/cube/")
 message('\n  Similarity analysis algorithm selected - Cube \n')
}else{
input.yaml$sim_dir <- as.character("sim_analyses/resnik/")
message('\n  Default to Resnik Similarity analysis algorithm \n')
}

yaml::write_yaml(input.yaml,paste0(input.yaml$output_dir,"input.yml") )

#Run sim_analysis files
if(input.yaml$sim_dir == as.character("sim_analyses/cube/")){
 source(paste0("scripts/cube_sim_auto_chunks.R"))
  source(paste0("gene_count_cube_auto.R"))
}
if(input.yaml$sim_dir == as.character("sim_analyses/resnik/")){
 source(paste0("scripts/res_mod_auto_chunks.R"))
  source(paste0("gene_count_resnik_mod_auto.R"))
} else{
  message('\n Please determine directory in input config file - Cant Proceed without that \n')
  break;
}
