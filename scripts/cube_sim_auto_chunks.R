library(tidyverse)

start <- Sys.time()
message(" \n Begin Cube chunking and comparison of patients \n ")


pat_table_prop <-  read_csv(paste0(input.yaml$output_dir, "pos_prop.csv"))

system(paste0("mkdir -p ",input.yaml$output_dir,input.yaml$sim_dir))
system(paste0("cp scripts/cube_master_sim_source.R ",input.yaml$output_dir,input.yaml$sim_dir,"."))
system(paste0("cp scripts/gene_count_cube_auto.R ",input.yaml$output_dir,input.yaml$sim_dir,"."))

setwd(paste0(input.yaml$output_dir,input.yaml$sim_dir))

#HPO terms in cohort (base)

#local_IC - use most recent version 
# local_IC <- read_csv("../../pos_IC.csv")
local_IC <- input.yaml$pos_ic %>% as.data.frame()

prop_temp <- pat_table_prop %>% 
  left_join(local_IC %>% select(HPO, prop.IC))

#################
#STEP 1: Set-up Intervals
#################
# subs = the amount of intervals you would like to split the patients into (height-wise)
##       the number of intervals = the number of parallel jobs that will be submitted into the cluster


pats = unique(pat_table_prop$famID)
n_pats = length(pats)
subs = input.yaml$n_subs# How many intervals you would like


#Rudimentary way to split the jobs into more equal sizes
v1 = c(5)
gap = 5

for(n in 2:subs){
  v1[n] = gap + floor((n_pats/subs)*(n-1))
}

v1[42] = n_pats
inter_subs = v1


#################
#STEP 2: Split cohort by intervals and submit to cluster
#################

#Split up the cohort (rowwise, to be submitted into the cluster in parallel)
  tot_sim_score=list()
  nstart = 1 

  for(n in 1:(length(inter_subs))){

    nend = inter_subs[n]

    destfile = paste("scn2a_cube_sub",n,".csv", sep = "")
    tot_sim_score[n] = destfile

    rname = paste("cube_sim_sub",n,".R", sep = "")

    simR = paste("sim_score = Compare_Cohort(prop_temp,", nstart, ",", nend,"); write.csv(sim_score," , "'",destfile,"',","row.names = F)",sep="")
    simsubR = paste("source('cube_master_sim_source.R');",simR) #Pasting the sourced R script which calculates the sim score
    write(simsubR,rname) #Creating a new R file

    shfile = paste0("/cm/shared/apps_chop/R/3.6.0/lib64/R/bin/Rscript ", rname)
    shname = paste("cube_sim_",n,".sh",sep = "")
    write(shfile,shname)
    system(paste0("qsub -cwd -l mem_free=5g,h_vmem=5g -V ",shname)) #submitting to cluster


    nstart = nend + 1

  }
  Sys.sleep(550)

#################
#STEP 3: Combine the sim_score intervals into one dataframe
#################

#Initalize comb_sim_scores
# 
destfile = paste0(tot_sim_score[1],"" )

  while (!file.exists(destfile)) {
    Sys.sleep(500)
  }
  
  if(file.exists(destfile)){
    f_name = "sub1"
    assign(f_name, read.csv(destfile, header=T),envir=globalenv())
    comb_sim_scores <- get(f_name)
  }


for(i in 2:length(tot_sim_score)){
  destfile = paste0(tot_sim_score[i],"" )

  while (!file.exists(destfile)) {
    Sys.sleep(240)
  }
  if(file.exists(destfile)){
    v_name = paste("sub", i, sep = "")
    assign(v_name, read.csv(destfile, header=T),envir=globalenv())
    n_scores <- get(v_name)
    comb_sim_scores <- rbind(comb_sim_scores,n_scores)
  }else{
    f_name = "sub1"
    assign(f_name, read.csv(destfile, header=T),envir=globalenv())
    comb_sim_scores <- get(f_name)
  }
}

comb_sim_scores[is.na(comb_sim_scores)] = 0
t_comb_sim_scores = t(comb_sim_scores)
diag(t_comb_sim_scores) = 0
t_comb_sim_scores[is.na(t_comb_sim_scores)] = 0

sim_score = t_comb_sim_scores + comb_sim_scores
diag(sim_score) = 0
# 
# 
write.csv(sim_score,"cube_sim_scn2a.csv", row.names = T) 
# 
#
# setwd(input.yaml$output_dir) 
message("\n  The Entire Cube chunking script ran successfully. Please find the output file mod_sim_scn2a.csv in the sim_analyses cube directory for further calculations of statistical significance \n ")
stop = Sys.time()
stop - start

