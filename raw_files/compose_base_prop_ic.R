#Compose Base and Prop

library(readr)
# library(readxl)
library(tidyverse)
library(Hmisc)

setwd("/Volumes/helbig_lab/projects/SCN2A/v15/")
# setwd("/mnt/isilon/helbig_lab/projects/SCN2A/v15/")

hpo_ancs <- read_csv("HPO_ancestors_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv") %>% 
  rename(HPO = term)

hpo_child <- read_csv("HPO_child_paths_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv") %>% 
  rename(HPO = term) %>% 
  mutate(children = paste0(HPO,";", children)) %>% 
  select(-X1)

hpo_neg_child <- hpo_child %>% 
  mutate(HPO = gsub("H","N",HPO)) %>% 
  mutate(children = gsub("H","N", children))
  

tot_base <- read_csv("SCN2A_full_V15.csv") %>% 
  # rename(famID = ID, HPO = term_id) %>% 
  # mutate(famID = paste0("fam",famID)) %>% 
  filter(variant_type_1 != "exclude") %>% 
  filter(grepl("HP:[[:digit:]]", HPO) | grepl("NP:[[:digit:]]", HPO) ) %>% 
  select(famID, HPO) %>% 
  unique()



write_csv(tot_base, "full_base_v15.csv")

########
# Base Files
########

neg_base <- tot_base %>% 
  filter(grepl("NP:[[:digit:]]", HPO)) %>% 
  select(famID, HPO) %>% 
  unique

pos_base <- tot_base %>% 
  filter(grepl("HP:[[:digit:]]", HPO)) %>% 
  unique

write_csv(neg_base, "neg_base_v15.csv")
write_csv(pos_base, "pos_base_v15.csv")

########
# Prop Files
########

# Remove modifier terms
redundant_terms <- c("HP:0003674","HP:0031797", "HP:0012823")

# "HP:0012823" ('clinical modifier'); "HP:0003674" ('Onset'); "HP:0031797" ('Clinical course')

pos_prop <- pos_base %>% 
  left_join(hpo_ancs) %>% 
  select(famID, ancs) %>% 
  separate_rows(ancs, sep=";") %>% 
  rename(HPO = ancs) %>% 
  filter(!is.na(HPO)) %>% 
  filter(HPO %nin% redundant_terms) %>%
  unique

neg_prop <- neg_base %>% 
  left_join(hpo_neg_child) %>% 
  select(famID, children) %>% 
  separate_rows(children, sep=";") %>% 
  rename(HPO = children) %>% 
  filter(!is.na(HPO), HPO !="NA") %>% 
  unique

write_csv(neg_prop, "neg_prop_v15.csv")
write_csv(pos_prop, "pos_prop_v15.csv")


full_prop <- neg_prop %>% rbind(pos_prop)

write_csv(full_prop, "full_prop_v15.csv")



#Negative Prop Pruned

neg_filter <- read_csv("/Volumes/helbig_lab/projects/SCN2A/v13/pruned_neg_filter_terms_v13.csv")

neg_prop_filt <- neg_prop %>% 
  filter(HPO %in% neg_filter$term)

write_csv(neg_prop_filt, "neg_prop_pruned_v15.csv")

########################
# Create IC
########################

########################################################################
########################################################################
########################################################################
########################################################################


neg_base_ct <- neg_base %>% 
  count(HPO) %>% 
  mutate(freq = n/length(unique(tot_base$famID)))

# pos_base = read_csv("/Volumes/helbig_lab/Dropbox/galerp/SCN2A/KC_SCN2A/scn2a_v2/scn2a_full_V4.csv") %>% 
#   rename(famID = ID, HPO = term_id) %>% 
#   # filter(!grepl("NP", HPO), grepl("HP:[[:digit:]]", HPO)) %>% 
#   filter(grepl("HP:[[:digit:]]", HPO)) %>% 
#   mutate(famID = paste0("fam_",famID)) %>% 
#   select(famID, HPO) %>% 
#   unique

pos_base_ct <- pos_base %>% 
  count(HPO) %>% 
  mutate(freq_base = n/length(unique(tot_base$famID))) %>% 
  mutate(base.IC = -log2(freq_base)) %>% 
  select(-n)

# hpo_ancs <- read_csv("/Volumes/helbig_lab/Dropbox/galerp/HPO_toolbox/HPO_ancestors_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv") %>% 
#   rename(HPO = term) 
# 
# hpo_child <- read_csv("/Volumes/helbig_lab/Dropbox/galerp/HPO_toolbox/HPO_child_paths_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv") %>% 
#   rename(HPO = term) %>% 
#   mutate(children = paste0(HPO,";",children)) %>% 
#   mutate(HPO = gsub("H","N", HPO)) %>% 
#   mutate(children = gsub("H","N", children)) %>% 
#   select(-X1)

pos_prop_ct <- pos_prop %>% 
  count(HPO) %>% 
  mutate(freq_prop = n/length(unique(tot_base$famID))) %>% 
  mutate(prop.IC = -log2(freq_prop)) %>% 
  select(-n)

pos_IC <- pos_base_ct %>% 
  full_join(pos_prop_ct) %>%
  right_join(hpo_ancs %>% select(HPO, def) %>% unique)

pos_IC[is.na(pos_IC)] <- 0

# tot_prop = pos_prop %>% rbind(neg_prop) %>% unique

length(unique(pos_prop$HPO))

write_csv(pos_IC,"pos_IC_v15.csv")


