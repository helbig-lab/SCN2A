#Using default method

#MICA algorithm
library(memoise)
library(tidyverse)

start <- Sys.time()
message(" \n Resnik mod chunking source file running... \n \n ")

setwd(input.yaml$sim_dir)

#Add IC Data

#local_IC - use most recent version 
local_IC <- read_csv(input.yaml$pos_ic) 


#HPO terms in cohort (base and propagated)
pat_table_base <- read_csv(paste0(input.yaml$file_path, "pos_base.csv"))
pat_table_prop <- read_csv(paste0(input.yaml$file_path, "pos_prop.csv"))


ancest_comp <- read_csv(paste0(input.yaml$file_path, "HPO_ancestors_v0.1.2.csv"))
path <- ancest_comp %>% 
  tidyr::separate_rows(ancs, sep = ";") %>% 
  filter(ancs != "") %>% 
  mutate(Counter1 = ancs) %>% 
  select(-def) %>% 
  rename(Term = term, HPO1 = ancs)
  



################ 
# MICA function
################ 
mica <- function(hpo1, hpo2){
  
  path1_unique <- path %>% dplyr::filter(Term == hpo1)
  # path1_unique <- path[which(path$Term == hpo1),]
  
  path2_unique <- path %>% dplyr::filter(Term == hpo2)
  # path2_unique <- path[which(path$Term == hpo2),]
  
  joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1')
  
  joint2 <- joint1 %>% left_join(local_IC, by = c('Counter1' = 'HPO')) 
  
  mica_ic <- joint2$prop.IC %>% max  
  
  
  return(mica_ic)
  
}

memo_mica <- memoise(mica)



################ 
#Patient Comparison Function
################ 

pat_compare <- function(pat1, pat2)
{
  hpo_pat1 <- pat_table_base %>% dplyr::filter(famID == pat1)
  hpo_pat2 <- pat_table_base %>% dplyr::filter(famID == pat2)
  
  # hpo_pat1 <- pat_table_base[which(pat_table_base$famID == pat1),]
  # hpo_pat2 <- pat_table_base[which(pat_table_base$famID == pat2),]
  
  
  #create data frame with HPO of pat1 in x, HPO of pat2 in y
  
  x_length <- length(hpo_pat1$HPO)
  
  y_length <- length(hpo_pat2$HPO)
  
  
  ic_matrix <- as.data.frame(matrix(ncol=x_length, nrow=y_length))
  
  names(ic_matrix) <- hpo_pat1$HPO
  
  rownames(ic_matrix) <- hpo_pat2$HPO
  
  for (i in 1:y_length){
    
    for(j in 1:x_length)
    {
      # HP = mixedsort(hpo_pat2$HPO[i],hpo_pat1$HPO[j])  #this will make the comparison memoise faster
      ic_matrix[i,j] <- memo_mica(hpo_pat2$HPO[i], hpo_pat1$HPO[j])
      
    }
  }
  
  
  max_col <- apply(ic_matrix,2,max)
  
  max_row <- apply(ic_matrix,1,max)
  
  max_complete <- sum(max_col,max_row)/2
  
  return(max_complete)
}


############### run on entire cohort 
#cohort_file is a prop file
Compare_Cohort=function(cohort_file,n1,n2){
  
  patients=unique(cohort_file$famID)
  
  dimension <- length(unique(patients))
  
  subcohort <- n2 - n1 + 1
  
  pat_matrix <- as.data.frame(matrix(ncol=dimension, nrow=subcohort))
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[n1:n2]
  for(x in 1:subcohort) #looping through each row
  {
    for(y in (x + n1 - 1):dimension) #looping through each column
    {
      # pat_matrix = sapply(1:dimension, function(y) (sapply(1:subcohort, function(x) pat_matrix[y,x] <- pat_compare(names(pat_matrix)[y], rownames(pat_matrix)[x]))))
      pat_matrix[x,y] <- pat_compare(names(pat_matrix)[y], rownames(pat_matrix)[x])
    }
  }
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[n1:n2]
  
  return(pat_matrix)
}

message("\n  The Resnik mod source script ran successfully. \n ")
stop = Sys.time()

