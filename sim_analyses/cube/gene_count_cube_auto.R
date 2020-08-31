######
#Load required libraries
######
library(tidyverse)
library(stringr)

start <- Sys.time()
message(" \n Begin Cube similarity analysis - test of statistical significance  \n ")

setwd(input.yaml$sim_dir)

'%!in%' <- function(x,y)!('%in%'(x,y))

print("start")
print(Sys.time())
######
#STEP 1: Load required files
######
#SCN2A dataset

paper_vars <- read_csv(paste0(input.yaml$file_path,"SCN2A_full.csv")) %>% 
  # rename(famID = ID, HPO = term_id) %>% 
  # mutate(famID = paste0("fam",famID)) %>% 
  filter(variant_type_1 != "exclude")

# seg_s4s5 <- paper_vars %>% 
#   select(famID, segment) %>% 
#   filter(segment == 'S4-S5') %>% 
#   unique() %>% 
#   mutate(var = 'S4-S5') %>% 
#   select(famID, var)

miss <- paper_vars %>% 
  filter(variant_type_2 == 'missense') %>% 
  mutate(var = "Missense") %>% 
  select(famID, var) %>% 
  unique

miss1 <- paper_vars %>% 
  filter(variant_type_2 == 'missense' & segment == 'S5-S6') %>% 
  mutate(var = 'Miss_S5-S6') %>% 
  select(famID, var) %>% 
  unique

miss2 <- paper_vars %>% 
  filter(variant_type_2 == 'missense' & segment != 'S5-S6') %>% 
  mutate(var = 'Miss_Not_S5-S6') %>% 
  select(famID, var) %>% 
  unique

ptv <- paper_vars %>% 
  filter(variant_type_2 == 'PTV') %>% 
  rename(var = variant_type_2) %>% 
  select(famID, var) %>% 
  unique

#Paper sources - Double check with Peter
paper <- paper_vars %>% 
  filter(!is.na(paper_source)) %>% 
  rename(var = paper_source) %>% 
  mutate(var = case_when(
    var == 'wolff 2017' ~ 'Wolff 2017',
    TRUE ~ var
  )) %>% select(famID, var) %>% 
  mutate(var = word(var, 1,2)) %>% 
  unique()

# R853Q <- paper_vars %>% 
#   select(famID, variant) %>% 
#   filter(variant == "R853Q") %>% 
#   rename(var = variant) %>% 
#   unique()

ptv_miss1 <- ptv %>% 
  rbind(miss1) %>% 
  mutate(var = 'S5-S6_and_PTV')

#All recurrent variants
rec_var <- paper_vars %>% 
  select(famID, variant) %>% 
  unique %>% 
  add_count(variant, sort=T) %>% 
  rename(var = variant) %>%
  filter(n > 1) %>% 
  select(famID, var) %>% 
  unique

#All domains
dom <- paper_vars %>% 
  filter(!is.na(domain)) %>% 
  filter(variant_type_2 == 'missense') %>% 
  select(famID, domain) %>% 
  unique() %>% 
  rename(var = domain) 

#All segments
seg <- paper_vars %>% 
  filter(!is.na(segment)) %>% 
  filter(variant_type_2 == 'missense') %>% 
  select(famID, segment) %>% 
  unique() %>% 
  rename(var = segment) %>% 
  select(famID, var)

#All variant classes

#All phenotypic groups
pheno_sub <- paper_vars %>%
  rename(var = broad_phx) %>%
  select(famID, var) %>%
  unique()

#GOF, LOF
gof_lof <- paper_vars %>%
  filter(!is.na(function.)) %>% 
  rename(var = function.) %>%
  select(famID, var) %>%
  unique()


variants <- miss %>% rbind(miss1, miss2, ptv, paper, ptv_miss1, rec_var, dom, seg, pheno_sub, gof_lof) %>% 
  unique

sim_score = read_csv(paste0(input.yaml$sim_dir,"cube_sim_scn2a.csv")) #%>% select(-c(1)) #%>% select(-X1) 
  sim_score<- do.call(data.frame, lapply(sim_score, function(x) {
  replace(x, is.infinite(x) | is.na(x), 0)
}))

diag(sim_score) <- 0
cor_names = names(sim_score)
names(sim_score) = cor_names
rownames(sim_score) = cor_names

###########################################################

######
#STEP 2: Clean and organize the data
##Generate required tables and dataframes
######
#Unpack fill with famIDs and their corresponding variants

famIDs_var = variants %>%
  dplyr::rename(variant = var) %>%
  separate_rows(famID, sep = ";") %>% 
  unique()


famIDs_sim <- sim_score %>% rownames %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>%  dplyr::mutate(sim = "sim_score")

fam_combined <- famIDs_sim %>%  dplyr::left_join(famIDs_var) %>% 
  filter(!is.na(sim), !is.na(variant)) %>% #only famIDs in the sim matrix (should be all)
  select(famID, sim, variant) %>% 
  unique #remove instances where patient has same gene which results in duplications


#DATA CHECK
##Trios without sim_scores or without variants
no_sim <- as.vector(famIDs_var$famID[famIDs_var$famID %!in% rownames(sim_score)]) #variant but no sim_score
no_var <- as.vector(rownames(sim_score)[rownames(sim_score) %!in% famIDs_var$famID]) #sim_score but no variant


######
#STEP 3: Isolate variants
######
# Table of denovo variants

#list of all genes that occur in more than one patient 
all_genes <- fam_combined %>%  dplyr::count(variant) %>% 
  dplyr::rename(Freq = n) %>% 
   filter(Freq > 1)
  # filter(Freq < 25)


#Create all 100k files- sending each in as a seperate job in the cluster
pat_pairs = all_genes$Freq %>% unique %>% as.list

if(input.yaml$cluster == 'c') {
  for( i in 1:length(pat_pairs)) {
    n = pat_pairs[i] %>% as.numeric()
    destfile = paste("100k_",n,"cube_scn2a.csv", sep = "")

    if(!file.exists(paste0(input.yaml$file_path,"sim_analyses", "/",destfile))){
      print("Error")
      print(destfile)
      ndraw = paste("n100k = sim_pat_draw(sim_score,", n, ");
                    write.csv(n100k,","'",destfile,"',","row.names = F)",sep="")
      rname = paste("cube_100k_",n,"p.R",sep = "")
      p100kR = paste("source('pat_100k_draw_npat_cube.R');",ndraw)
      write(p100kR,rname)

      shfile = paste0(input.yaml$sim_cluster, " ", rname)
      shname = paste("p",n,"_100k_cube.sh",sep = "")
      write(shfile,shname)
      system(paste0("qsub -cwd -l mem_free=4g,h_vmem=4g -V ",shname))
    }
    }
} else {
  for( i in 1:length(pat_pairs)) {
    n = pat_pairs[i] %>% as.numeric()
    destfile = paste("100k_",n,"cube_scn2a.csv", sep = "")

    if(!file.exists(paste0(input.yaml$file_path,"sim_analyses", "/",destfile))){
      print("Error")
      print(destfile)
      ndraw = paste("n100k = sim_pat_draw(sim_score,", n, ");
                    write.csv(n100k,","'",destfile,"',","row.names = F)",sep="")
      rname = paste("cube_100k_",n,"p.R",sep = "")
      p100kR = paste("source('pat_100k_draw_npat_cube.R');",ndraw)
      write(p100kR,rname)

      shfile = paste0(input.yaml$sim_cluster, " ", rname)
      shname = paste("p",n,"_100k_cube.sh",sep = "")
      write(shfile,shname)
      system(shname)
    }
    }
}


########
#Function: gene_df
########
#Input: gene of interest from the data set
#Output: return table of famIDs with gene compared to one another (i.e. with their sim scores) 

gene_df <- function(gene)
{
  az <- fam_combined %>% filter(variant==gene) %>% unique
  len1 = nrow(az)
  
  #matrix of all combinations
  matrix = t(combn(1:len1, 2))
  len2 = nrow(matrix)
  
  df <- as.data.frame(matrix(ncol=4,nrow=len2))
  names(df) <- c('fam1','fam2','gene','sim_score')
  df[,3] <- gene
  
  #fill df with indices of matrix and add sim_score
  for (i in 1:len2)
  {
    x_name = az$famID[matrix[i,1]]
    y_name = az$famID[matrix[i,2]]
    df[i,1] <- x_name
    df[i,2] <- y_name
    #y_name_cor <- gsub("-",".",y_name) #column names in dataframe replace "_" with "."
    df[i,4] <- sim_score[x_name,y_name]
    if(is.na(df[i,4])){
      print(x_name)
      print(y_name)
    }
  }
  df <- df[!duplicated(df),] #matrix of unique combinations
  return(df)
}

#Creating dataframe of similarity comparisons with every combination of patient pairs 
##with the same denove gene

#initialize with first gene
pair_corrected <- gene_df(all_genes$variant[1])

#Create for all genes
for (i in 2:nrow(all_genes)){
print(all_genes$variant[i])
  temp <- gene_df(all_genes$variant[i])
  pair_corrected <- rbind(pair_corrected,temp)
}
pair_corrected <- pair_corrected %>% filter(!is.na(sim_score))

#########
#STEP 4: Create similarity score distributions based on patients with the same gene
## Creating large dataframe of pvalues based on similarity between patients with
###the same denovo gene
#######

#Function: sim_random
##Output: picks sim_score from two random pairs
sim_random <- function() {
  random = sample(1:nrow(sim_score),2,replace = FALSE)
  sims_rnd <- sim_score[random[1],random[2]]
  return(sims_rnd)
}

#Function: draw_av
##Output: returns the average sim_score in n random individuals 
draw_av <- function(n) {
  count = 0
  for (i in 1:n) {
    count = count + sim_random()
  }
  count = count/n
  return(count)
}

#Function: draw_median
##Output: returns the median sim_score in n random individuals 
draw_median <- function(n) {
  count = c()
  for (i in 1:n) {
    count[i] = sim_random()
  }
  median_i <- median(count)
  return(median_i)
}

#Function: draw_median
##Output: returns the mode sim_score in n random individuals 
estimate_mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}

#Function: draw_mode
##Output: returns the mode sim_score in n random individuals 
draw_mode <- function(n) {
  count = c()
  for (i in 1:n) {
    count[i] = sim_random()
  }
  mode_i <- estimate_mode(count)
  return(mode_i)
}


#Determine number of pairs per gene
gene_x <- unique(pair_corrected$gene)

gene_count <- as.data.frame(matrix(ncol=9,nrow=length(gene_x))) 

names(gene_count) <- c("gene","n_pats","pairs","av_sim","median_sim","mode_sim","p_av","p_median","p_mode")

gene_count <- gene_count %>% mutate(gene = gene_x)

#Finding the average, median, and mode similarity,for each gene in the cohort
for (i in 1:nrow(gene_count)) {
  name_x <- subset(pair_corrected,pair_corrected$gene == as.character(gene_count[i,c("gene")]))
  gene_count[i,c('n_pats')] <- vctrs::vec_c(name_x$fam1, name_x$fam2) %>% unique %>% length()
  gene_count[i,c('pairs')] <- nrow(name_x)  
  gene_count[i,c('av_sim')] <- sum(name_x$sim_score)/nrow(name_x)
  gene_count[i,c('median_sim')] <- median(name_x$sim_score)
  if (length(name_x$sim_score) > 1) {gene_count[i,c('mode_sim')] <- estimate_mode(name_x$'sim_score')}
  if (length(name_x$sim_score) == 1) {gene_count[i,c('mode_sim')] <- name_x$sim_score}
}

#########
#STEP 5: Find the p_value for each gene's similarity
## using the average, median, and mode similarity,for each gene from gene_count table
#########


#Function: exact_p_av2
##Determines p-value for the average
exact_p_av2 <- function(sample_size,sim_score_med,number_random_draws)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  v_name <- paste("n",n,"_","100k",sep = "")
  
  if(exists(v_name)){
    ax <-get(v_name)
    
  }else  {
    destfile = paste("100k_",n,"cube_scn2a.csv", sep = "")
    if(file.exists(destfile)){
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }else{
      while (!file.exists(destfile)) {
        Sys.sleep(300)
      }
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }
  }
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x){

    foo = ecdf(ax$mean)
    return(1 - foo(x))
    
  }
  
  p <- p_value(sim_score_med)
  return(p)
}


#Function: exact_p_med2
##Determines p-value for the median
exact_p_med2 <- function(sample_size,sim_score_med,number_random_draws)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  v_name <- paste("n",n,"_","100k",sep = "")
  
  if(exists(v_name)){
    ax <-get(v_name)
    
  }else  {
    destfile = paste("100k_",n,"cube_scn2a.csv", sep = "")
    if(file.exists(destfile)){
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }else{
      while (!file.exists(destfile)) {
        Sys.sleep(300)
      }
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }
  }
  
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {

    foo = ecdf(ax$median)
    return(1 - foo(x))
    
  }
  p <- p_value(sim_score_med)
  return(p)
  
}

#Function: exact_p_mod2
##Determines p-value for the mode
exact_p_mod2 <- function(sample_size,sim_score_mod,number_random_draws){
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  v_name <- paste("n",n,"_","100k",sep = "")
  
  if(exists(v_name)){
    ax <-get(v_name)
    
  }else  {
    destfile = paste("100k_",n,"cube_scn2a.csv", sep = "")
    if(file.exists(destfile)){
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }else{
      while (!file.exists(destfile)) {
        Sys.sleep(300)
      }
      assign(v_name, read.csv(destfile, header=T),envir=globalenv())
      ax <- get(v_name)
    }
  }

  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {

    foo = ecdf(ax$mode)
    return(1 - foo(x))
    
  }
  p <- p_value(sim_score_mod)
  return(p)
}

#loop through each gene for p_av, p_med
for (i in 1:nrow(gene_count)){

  e1 = gene_count$n_pats[i] #sample size
  #e2 = gene_count$mode_sim[i] #sim score mod
  e3 = 100000 # 10K unless n=1, then 100K (predefined)
  p_average = exact_p_av2(e1,gene_count$av_sim[i],e3)
  p_median = exact_p_med2(e1,gene_count$median_sim[i],e3)
  p_mod = exact_p_mod2(e1,gene_count$mode_sim[i],e3)
  gene_count[i,7] <- p_average
  gene_count[i,8] <- p_median
  gene_count[i,9] <- p_mod
  
}


write.csv(gene_count,paste0(input.yaml$output_dir,"sim_analyses/gene_count_cube_scn2a.csv"),row.names = F)

message("\n  Cube similarity analysis - test of statistical significance ran successfully. \n ")
stop = Sys.time()
stop-start
