library(tidyverse)

start <- Sys.time()
message(" \n Begin Resnik chunking and comparison of patients \n ")

# setwd(input.yaml$sim_dir)

#HPO terms in cohort (base)

pat_table_base <- read_csv(paste0(input.yaml$output_dir, "pos_base.csv"))
pat_table_prop <- read_csv(paste0(input.yaml$output_dir, "pos_prop.csv"))

system(paste0("mkdir -p ",input.yaml$output_dir,input.yaml$sim_dir))
system(paste0("cp scripts/sim_res_mod_auto_source.R ",input.yaml$output_dir,input.yaml$sim_dir,"."))
system(paste0("cp scripts/gene_count_resnik_mod_auto.R ",input.yaml$output_dir,input.yaml$sim_dir,"."))
setwd(paste0(input.yaml$output_dir,input.yaml$sim_dir))
#################
#STEP 1: Set-up Intervals
#################
# subs = the amount of intervals you would like to split the patients into (height-wise)
##       the number of intervals = the number of parallel jobs that will be submitted into the cluster

pats = unique(pat_table_prop$famID)
n_pats = length(pats)
subs = input.yaml$n_subs  # How many intervals you would like


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

    destfile = paste0("mod_scn2a_sub",n,".csv")
    tot_sim_score[n] = destfile

    rname = paste("mod_sub",n,".R", sep = "")

    simR = paste("sim_score = Compare_Cohort(pat_table_prop,", nstart, ",", nend,"); write.csv(sim_score," , "'",destfile,"',","row.names = F)",sep="")
    simsubR = paste0("source('sim_res_mod_auto_source.R');",simR) #Pasting the sourced R script which calculates the sim score
    write(simsubR,rname) #Creating a new R file

    shfile = paste0("/cm/shared/apps_chop/R/3.6.0/lib64/R/bin/Rscript ", rname)
    shname = paste("mod_",n,".sh",sep = "")
    write(shfile,shname)
    system(paste0("qsub -cwd -l mem_free=5g,h_vmem=5g -V ",shname)) #submitting to cluster


    nstart = nend + 1

  }


#################
#STEP 3: Combine the sim_score intervals into one dataframe
#################
Sys.sleep(900) #wait 15 minutes
#Initalize comb_sim_scores
# 
destfile = paste0(tot_sim_score[1],"" )
if(file.exists(destfile)){ #check if first sim sub file created
  f_name = "sub1"
  assign(f_name, read.csv(destfile, header=T),envir=globalenv())
  comb_sim_scores <- get(f_name)

}else{
  while (!file.exists(destfile)) {
    Sys.sleep(450)
  }
  f_name = "sub1"
  assign(f_name, read.csv(destfile, header=T),envir=globalenv())
  comb_sim_scores <- get(f_name)
}

for(i in 2:length(tot_sim_score)){
  destfile = paste0(tot_sim_score[i],"" )


  if(file.exists(destfile)){
    v_name = paste("sub", i, sep = "")
    assign(v_name, read.csv(destfile, header=T),envir=globalenv())
    n_scores <- get(v_name)
    comb_sim_scores <- rbind(comb_sim_scores,n_scores)
  }else{
    while (!file.exists(destfile)) {
      Sys.sleep(240)
    }
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
# 
# 
write_csv(sim_score,paste0("mod_sim_scn2a.csv"))
# 
# 
# setwd(input.yaml$output_dir)
message("\n  The Entire Resnik chunking script ran successfully. Please find the output file mod_sim_scn2a.csv in the sim_analyses resnik directory for further calculations of statistical significance \n ")
stop = Sys.time()
message(paste(stop - start))

