start <- Sys.time()
message(" \n Starting config file... \n ")


message("\n  ...Checking for analyses to run... \n ")

#Term propagation
if(is.null(input.yaml$pos_ic) == F ){
  next;
}
else{
    message("\n  Manually propagating information content... \n ")
    source("raw_files/compose_base_prop_ic.R))
    }
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


