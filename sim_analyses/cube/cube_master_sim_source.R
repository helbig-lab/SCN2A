library(tidyverse)

start <- Sys.time()
message(" \n Cube source file running... \n ")

pat_table_prop <- read_csv(paste0(input.yaml$file_path, "pos_prop_v14.csv"))

#local_IC - use most recent version 
local_IC <- read_csv(paste0(input.yaml$file_path,"pos_IC_v14.csv")) 

prop_temp <- pat_table_prop %>% 
  left_join(local_IC %>% select(HPO, prop.IC))


Compare_Cohort=function(prop_temp, n1, n2){
  
  
  patients = unique(prop_temp$famID)
  
  dimension <- length(unique(patients))
  
  subcohort <- n2 - n1 + 1
  
  pat_matrix <- as.data.frame(matrix(ncol=dimension, nrow=subcohort))
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[n1:n2] 
  for(x in 1:subcohort) #looping through each row
  {
    for(y in (x + n1 - 1):dimension) #looping through each column
    {
      
      pat1 <- prop_temp %>% filter(famID == names(pat_matrix)[y])  #patients[x])  
      pat2 <- prop_temp %>% filter(famID == rownames(pat_matrix)[x])     #patients[y]) 
      
      pat_comp <- pat1 %>% inner_join(pat2, by = c("HPO","prop.IC"))
      
      pat_matrix[x,y] <- sum(pat_comp$prop.IC)
    }
  }
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[n1:n2]
  
  return(pat_matrix)
}

message("\n  Cube source file completed \n ")
stop = Sys.time()
stop - start



