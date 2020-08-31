library(tidyverse)

start <- Sys.time()
message(" \n Begin Resnik iteration source file \n ")

setwd(input.yaml$sim_dir)


sim_score = read_csv(paste0(input.yaml$sim_dir,"mod_sim_scn2a.csv"))
sim_score<- do.call(data.frame, lapply(sim_score, function(x) {
  replace(x, is.infinite(x) | is.na(x), 0)
})
)

diag(sim_score) <- 0
cor_names = names(sim_score)
names(sim_score) = cor_names
rownames(sim_score) = cor_names



estimate_mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}


sim_pat_draw = function(sim_score, num_pats)  {
  
  r_100k = as.data.frame(matrix(nrow = 100000, ncol = 3))
  names(r_100k) = c("median","mean", "mode")
  
  pat_vect = names(sim_score)
  
  for(n in 1: nrow(r_100k)){
    IDs = sample(pat_vect, num_pats)
    sub_sim = sim_score[(rownames(sim_score) %in% IDs), (names(sim_score) %in% IDs)]
    diag(sub_sim) = 12345
    vect_scores =  unlist(sub_sim)
    vect_scores = vect_scores[-which(vect_scores == 12345)]
    
    r_100k$median[n] = median(vect_scores)
    r_100k$mean[n] = mean(vect_scores)
    r_100k$mode[n] = estimate_mode(vect_scores)
  }
  return(r_100k)
}

message("\n  The Resnik iteration source filescript ran successfully. \n ")
stop = Sys.time()


