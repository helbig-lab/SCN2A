library(tidyverse)

start <- Sys.time()
message(" \n Starting frequency anaylses... \n ")

######
# Read in files
######

variants <- read.csv(paste0(input.yaml$file_path, "SCN2A_full.csv"), 
                     stringsAsFactors = FALSE) %>% 
  filter(!is.na(famID)) %>% 
  filter(variant_type_1 != "exclude") %>% 
  filter(grepl("HP:[0-9]", HPO)) %>% 
  unique()

hpo_def <- read.csv(paste0(input.yaml$file_path, "HPO_ancestors_v0.1.2.csv"), 
                    stringsAsFactors = FALSE) %>% 
  rename(HPO = term) %>% 
  select(HPO, def) %>% 
  unique()

# Positive propagated terms
prop <- read.csv(paste0(input.yaml$file_path, "pos_base.csv"), 
                 stringsAsFactors = FALSE)

# Frequency - all
ic <- read.csv(paste0(input.yaml$file_path,"pos_IC.csv"), 
               stringsAsFactors = F) %>% select(HPO, freq_prop)


######
# Clean up and merge
######

# Variant class groups
variants$variant_type_3 <- case_when(variants$variant_type_2 == 'PTV' ~ 'PTV',
                                     variants$variant_type_2 == 'missense' & variants$segment == 'S5-S6' ~ 'missense S5-S6',
                                     variants$variant_type_2 == 'missense' & variants$segment != 'S5-S6' ~ 'other missense')

# Merge with prop terms
merged <- prop %>% left_join(variants %>% select(famID, broad_phx, 
                                                 variant_type_2, variant_type_3, 
                                                 domain, segment) %>% unique())
unique_terms <- merged$HPO %>% unique()

######
# Table 1 - frequencies only
######

## Phenotypic groups
pheno_groups <- variants$broad_phx %>% unique()
pheno_groups <- pheno_groups[c(2, 4, 1, 5, 3)]
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(pheno_groups)) {
  
  pg = pheno_groups[i]
  
  pg_terms <- merged %>% filter(broad_phx == pg)
  pats = pg_terms$famID %>% unique()
  
  term_freq <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(pg == "encephalopathy" ~ "DEE",
                                                          pg == "ASD" ~ "ASD",
                                                          pg == "benign" ~ "BFNIS",
                                                          pg == "epilepsy" ~ "Other_EPI",
                                                          pg == "atypical" ~ "ATYPICAL")

}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"phenotypic_groups_frequencies.csv"), row.names = FALSE)


## Variant classes - "variant_type_2"
variant_groups_2 <- variants %>% filter(variant_type_2 != "") %>% select(variant_type_2) %>% pull() %>% unique()
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_2)) {

  pg = variant_groups_2[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))

  pg_terms <- merged_1 %>% filter(variant_type_2 == pg)
  pats = pg_terms$famID %>% unique()

  term_freq <- pg_terms %>%
    dplyr::group_by(HPO) %>%
    dplyr::summarize(freq=n()/length(pats))

  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0

  names(freq_all)[names(freq_all) == 'freq'] <- case_when(pg == "missense" ~ "missense",
                                                          pg == "PTV" ~ "PTV")

}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"ptv_missense_frequencies.csv"), row.names = FALSE)


## Variant classes - "variant_type_3" (S5-S6)
variant_groups_3 <- variants$variant_type_3 %>% na.exclude %>% unique()
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_3)) {
  
  pg = variant_groups_3[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_3))
  
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  
  term_freq <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(pg == "missense S5-S6" ~ "missense S5-S6",
                                                          pg == "other missense" ~ "missense other",
                                                          pg == "PTV" ~ "PTV")
  
}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"ptv_missense_S5.S6_frequencies.csv"), row.names = FALSE)


## Variant classes - "variant_type_4" - Recurrent
variants %>% select(variant, famID) %>% unique() %>% group_by(variant) %>% count() %>% filter(n>1) -> recur

rec_variants <- variants %>% filter(variant %in% recur$variant) %>% select(famID, variant) %>% unique()

merged_2 <- prop %>% left_join(rec_variants %>% select(famID, variant) %>% unique())
unique_terms <- merged$HPO %>% unique()

recur_groups <- rec_variants$variant %>% unique()
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(recur_groups)) {
  
  rg = recur_groups[i]

  rg_terms <- merged_2 %>% filter(variant == rg)
  pats = rg_terms$famID %>% unique()
  
  term_freq <- rg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(TRUE ~ rg)
  
}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"recurrent_var_frequencies.csv"), row.names = FALSE)


## Variant classes - "Miss_Not_S5_S6"
variant_groups_4 <- variants$variant_type_3[variants$variant_type_3 == "other missense"] %>% na.exclude %>% unique()
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_4)) {
  
  pg = variant_groups_4[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_3))
  
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  
  term_freq <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(pg == "other missense" ~ "Miss_Not_S5_S6")
  
}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"missense_not_S5.S6_frequencies.csv"), row.names = FALSE)


## Segment groups
segment_groups <- variants %>% filter(variant_type_2 == 'missense') %>% 
  select(segment) %>% na.exclude %>% unique()
segment_groups <- segment_groups[-c(3)]
segment_groups <- as.vector(segment_groups[,1])
freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(segment_groups)) {
  
  sg = segment_groups[i]

  sg_terms <- merged %>% filter(segment == sg)
  #pats = sg_terms$famID %>% unique()
  pats = sg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  
  term_freq <- sg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(TRUE ~ sg)
  
}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"segment_frequencies.csv"), row.names = FALSE)


## domain groups
domain_groups <- variants %>% filter(variant_type_2 == 'missense') %>% 
  select(domain) %>% na.exclude %>% unique()
domain_groups <- domain_groups[-c(3)]
domain_groups <- as.vector(domain_groups[,1])

freq_all <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(domain_groups)) {
  
  dg = domain_groups[i]
  
  dg_terms <- merged %>% filter(domain == dg)
  #pats = dg_terms$famID %>% unique()
  pats = dg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  
  term_freq <- dg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- case_when(TRUE ~ dg)
  
}

freq_all <- freq_all %>% left_join(ic) %>% rename(ALL = freq_prop)
freq_all <- freq_all[order(-freq_all$ALL),]

write.csv(freq_all, paste0(input.yaml$freq_dir,"domain_frequencies.csv"), row.names = FALSE)


keep(variants, merged, ic, hpo_def, unique_terms, pheno_groups, variant_groups_2, 
     variant_groups_3, variant_groups_4, recur_groups, merged_2, segment_groups, sure = T)



######
# Table 2 - raw counts and frequency stats
######

## Phenotypic groups
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(pheno_groups)) {
  
  pg = pheno_groups[i]
  
  pg_terms <- merged %>% filter(broad_phx == pg)
  pats = pg_terms$famID %>% unique()
  
  term_count <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(pg == "encephalopathy" ~ "DEE_present",
                                                              pg == "ASD" ~ "ASD_present",
                                                              pg == "benign" ~ "BFNIS_present",
                                                              pg == "epilepsy" ~ "Other_EPI_present",
                                                              pg == "atypical" ~ "ATYPICAL_present")
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(pg == "encephalopathy" ~ "DEE_absent",
                                                              pg == "ASD" ~ "ASD_absent",
                                                              pg == "benign" ~ "BFNIS_absent",
                                                              pg == "epilepsy" ~ "Other_EPI_absent",
                                                              pg == "atypical" ~ "ATYPICAL_absent")

  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(pg == "encephalopathy" ~ "DEE_f",
                                                          pg == "ASD" ~ "ASD_f",
                                                          pg == "benign" ~ "BFNIS_f",
                                                          pg == "epilepsy" ~ "Other_EPI_f",
                                                          pg == "atypical" ~ "ATYPICAL_f")
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(pg == "encephalopathy" ~ "DEE_lower",
                                                              pg == "ASD" ~ "ASD_lower",
                                                              pg == "benign" ~ "BFNIS_lower",
                                                              pg == "epilepsy" ~ "Other_EPI_lower",
                                                              pg == "atypical" ~ "ATYPICAL_lower")
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(pg == "encephalopathy" ~ "DEE_upper",
                                                                 pg == "ASD" ~ "ASD_upper",
                                                                 pg == "benign" ~ "BFNIS_upper",
                                                                 pg == "epilepsy" ~ "Other_EPI_upper",
                                                                 pg == "atypical" ~ "ATYPICAL_upper")
  
  
}

raw_counts <- raw_counts[order(-raw_counts$DEE_f),]

write.csv(raw_counts, "frequency_analyses/phenotypic_groups_raw_counts.csv", row.names = FALSE)


## Variant classes - "variant_type_2"
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_2)) {
  
  pg = variant_groups_2[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))
  
  pg_terms <- merged_1 %>% filter(variant_type_2 == pg)
  pats = pg_terms$famID %>% unique()
  
  term_count <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(pg == "missense" ~ "missense_present",
                                                                   pg == "PTV" ~ "PTV_present")
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(pg == "missense" ~ "missense_absent",
                                                                  pg == "PTV" ~ "PTV_absent")
  
  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(pg == "missense" ~ "missense_f",
                                                              pg == "PTV" ~ "PTV_f")
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(pg == "missense" ~ "missense_lower",
                                                                 pg == "PTV" ~ "PTV_lower")
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(pg == "missense" ~ "missense_upper",
                                                                 pg == "PTV" ~ "PTV_upper")
  
  
}

raw_counts <- raw_counts[order(-raw_counts$missense_f),]

write.csv(raw_counts, "frequency_analyses/ptv_missense_raw_counts.csv", row.names = FALSE)


## Variant classes - "variant_type_3"
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_3)) {
  
  pg = variant_groups_3[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_3))
  
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  
  term_count <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(pg == "other missense" ~ "other_missense_present",
                                                                   pg == "missense S5-S6" ~ "S5.S6_missense_present",
                                                                   pg == "PTV" ~ "PTV_present")
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(pg == "other missense" ~ "other_missense_absent",
                                                                  pg == "missense S5-S6" ~ "S5.S6_missense_absent",
                                                                  pg == "PTV" ~ "PTV_absent")
  
  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(pg == "other missense" ~ "other_missense_f",
                                                              pg == "missense S5-S6" ~ "S5.S6_missense_f",
                                                              pg == "PTV" ~ "PTV_f")
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(pg == "other missense" ~ "other_missense_lower",
                                                                 pg == "missense S5-S6" ~ "S5.S6_missense_lower",
                                                                 pg == "PTV" ~ "PTV_lower")
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(pg == "other missense" ~ "other_missense_upper",
                                                                 pg == "missense S5-S6" ~ "S5.S6_missense_upper",
                                                                 pg == "PTV" ~ "PTV_upper")
  
  
}

raw_counts <- raw_counts[order(-raw_counts$other_missense_f),]

write.csv(raw_counts, "frequency_analyses/ptv_missense_S5.S6_raw_counts.csv", row.names = FALSE)


## Variant classes - "missense_not_S5_S6"
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(variant_groups_4)) {
  
  pg = variant_groups_4[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_3))
  
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  
  term_count <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(pg == "other missense" ~ "other_missense_present",
                                                                   pg == "missense S5-S6" ~ "S5.S6_missense_present",
                                                                   pg == "PTV" ~ "PTV_present")
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(pg == "other missense" ~ "other_missense_absent",
                                                                  pg == "missense S5-S6" ~ "S5.S6_missense_absent",
                                                                  pg == "PTV" ~ "PTV_absent")
  
  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(pg == "other missense" ~ "other_missense_f",
                                                              pg == "missense S5-S6" ~ "S5.S6_missense_f",
                                                              pg == "PTV" ~ "PTV_f")
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(pg == "other missense" ~ "other_missense_lower",
                                                                 pg == "missense S5-S6" ~ "S5.S6_missense_lower",
                                                                 pg == "PTV" ~ "PTV_lower")
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(pg == "other missense" ~ "other_missense_upper",
                                                                 pg == "missense S5-S6" ~ "S5.S6_missense_upper",
                                                                 pg == "PTV" ~ "PTV_upper")
  
  
}

raw_counts <- raw_counts[order(-raw_counts$other_missense_f),]

write.csv(raw_counts, "frequency_analyses/missense_not_s5_s6_raw_counts.csv", row.names = FALSE)


## Recurrent Variant groups
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(recur_groups)) {
  
  rg = recur_groups[i]
  
  rg_terms <- merged_2 %>% filter(variant == rg)
  pats = rg_terms$famID %>% unique()
  
  term_count <- rg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
 
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(TRUE ~ paste0(rg,"_present"))

  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(TRUE ~ paste0(rg,"_absent"))

  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(TRUE ~ paste0(rg,"_freq"))

  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(TRUE ~ paste0(rg,"_f_lower"))

  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(TRUE ~ paste0(rg,"_f_upper"))
  
  
}

write.csv(raw_counts, "frequency_analyses/recurrent_var_raw_counts.csv", row.names = FALSE)


## Segment groups
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(segment_groups)) {
  
  sg = segment_groups[i]
  
  sg_terms <- merged %>% filter(variant_type_2 == 'missense') %>% filter(segment == sg)
  #pats = sg_terms$famID %>% unique()
  pats = sg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  
  term_count <- sg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(TRUE ~ paste0(sg,"_present"))
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(TRUE ~ paste0(sg,"_absent"))
  
  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(TRUE ~ paste0(sg,"_f"))
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(TRUE ~ paste0(sg, "_lower"))
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(TRUE ~ paste0(sg,"_upper"))
  
  
}

write.csv(raw_counts, "frequency_analyses/segment_groups_raw_counts.csv", row.names = FALSE)

## Domain groups
raw_counts <- hpo_def %>% filter(HPO %in% unique_terms)

for (i in 1:length(domain_groups)) {
  
  dg = domain_groups[i]
  
  dg_terms <- merged %>% filter(variant_type_2 == 'missense') %>% filter(domain == dg)
  #pats = dg_terms$famID %>% unique()
  pats = dg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  
  term_count <- dg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(n_present = n(), 
                     n_absent = length(pats) - n(), 
                     freq = n()/length(pats))
  
  # Prop test for 95% CI for freq
  term_count$f_lower = NA
  term_count$f_upper = NA
  for (j in 1:nrow(term_count)) {
    prop_test = prop.test(term_count$n_present[j], (term_count$n_present[j] + term_count$n_absent[j]))
    term_count$f_lower[j] <- prop_test$conf.int[1]
    term_count$f_upper[j] <- prop_test$conf.int[2]
  }
  
  raw_counts <- raw_counts %>% left_join(term_count)
  raw_counts[is.na(raw_counts)] = 0
  
  names(raw_counts)[names(raw_counts) == 'n_present'] <- case_when(TRUE ~ paste0(dg,"_present"))
  
  names(raw_counts)[names(raw_counts) == 'n_absent'] <- case_when(TRUE ~ paste0(dg,"_absent"))
  
  names(raw_counts)[names(raw_counts) == 'freq'] <- case_when(TRUE ~ paste0(dg,"_f"))
  
  names(raw_counts)[names(raw_counts) == 'f_lower'] <- case_when(TRUE ~ paste0(dg, "_lower"))
  
  names(raw_counts)[names(raw_counts) == 'f_upper'] <- case_when(TRUE ~ paste0(dg,"_upper"))
  
  
}

write.csv(raw_counts, "frequency_analyses/domain_groups_raw_counts.csv", row.names = FALSE)




keep(variants, merged, ic, hpo_def, unique_terms, pheno_groups, variant_groups_2, variant_groups_3, recur_groups, merged_2, sure = T)

######
# Table 3 - Fisher's test phenotypic groups
######

sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','PHENO_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the phenotypic groups
for (i in 1:length(pheno_groups)) {
  
  pg = pheno_groups[i]
  
  # In phenotypic group
  pg_terms <- merged %>% filter(broad_phx == pg)
  pats = pg_terms$famID %>% unique()
  yes_hpo = merged %>% filter(famID %in% pats)
  
  # Not in phenotypic group
  no_hpo <- merged %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','PHENO_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$PHENO_GROUP <- pg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)
sig_test$PHENO_GROUP <- case_when(sig_test$PHENO_GROUP == "encephalopathy" ~ "DEE",
                                  sig_test$PHENO_GROUP == "ASD" ~ "ASD",
                                  sig_test$PHENO_GROUP == "benign" ~ "BFNIS",
                                  sig_test$PHENO_GROUP == "epilepsy" ~ "Other_EPI",
                                  sig_test$PHENO_GROUP == "atypical" ~ "ATYPICAL")

## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/phenotypic_groups_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/phenotypic_groups_sig_assoc_POS_OR.csv", row.names = FALSE)

  

######
# Table 3 - Fisher's test variant missense vs PTV groups
######


sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the variant groups
for (i in 1:length(variant_groups_2)) {
  
  pg = variant_groups_2[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))
  
  # In phenotypic group
  pg_terms <- merged_1 %>% filter(variant_type_2 == pg)
  pats = pg_terms$famID %>% unique()
  yes_hpo = merged_1 %>% filter(famID %in% pats)
  
  # Not in phenotypic group
  no_hpo <- merged_1 %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$VAR_GROUP <- pg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)


## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/ptv_missense_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/ptv_missense_sig_assoc_POS_OR.csv", row.names = FALSE)

######
# Table 3 - Fisher's test variant missense vs PTV groups WITH S5-S6
######


sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the variant groups
for (i in 1:length(variant_groups_3)) {
  
  pg = variant_groups_3[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))
  
  # In phenotypic group
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  yes_hpo = merged_1 %>% filter(famID %in% pats)
  
  # Not in phenotypic group
  no_hpo <- merged_1 %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$VAR_GROUP <- pg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)


## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/ptv_missense_S5.S6_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/ptv_missense_S5.S6_sig_assoc_POS_OR.csv", row.names = FALSE)


######
# Table 3 - Fisher's test variant missense NOT S5-S6
######


sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the variant groups
for (i in 1:length(variant_groups_4)) {
  
  pg = variant_groups_4[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_3))
  
  # In phenotypic group
  pg_terms <- merged_1 %>% filter(variant_type_3 == pg)
  pats = pg_terms$famID %>% unique()
  yes_hpo = merged_1 %>% filter(famID %in% pats)
  
  # Not in phenotypic group
  no_hpo <- merged_1 %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$VAR_GROUP <- pg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)


## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/missense_not_S5.S6_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/missense_not_S5.S6_sig_assoc_POS_OR.csv", row.names = FALSE)


############
#Table 3 - Recurrent Variants
############

sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','SEGMENT_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the recurrent groups
for (i in 1:length(recur_groups)) {
  
  rg = recur_groups[i]
  
  # In recurrent group
  rg_terms <- merged_2 %>% filter(variant == rg)
  pats = rg_terms$famID %>% unique()
  yes_hpo = merged_2 %>% filter(famID %in% pats)
  
  # Not in recurrent group
  no_hpo <- merged %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','RECCUR_VAR', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$RECCUR_VAR <- rg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)

## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/recurrent_var_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/recurrent_var_sig_assoc_POS_OR.csv", row.names = FALSE)


######
# Table 3 - Fisher's test segment groups
######

sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','SEGMENT_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the phenotypic groups
for (i in 1:length(segment_groups)) {
  
  sg = segment_groups[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))
  
  # In segment group
  sg_terms <- merged_1 %>% filter(segment == sg)
  #pats = sg_terms$famID %>% unique()
  pats = sg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  yes_hpo = merged_1 %>% filter(famID %in% pats)
  
  # Not in segment group
  no_hpo <- merged_1 %>% filter(variant_type_2 == 'missense') %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','SEGMENT_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$SEGMENT_GROUP <- sg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)
sig_test$SEGMENT_GROUP <- case_when(TRUE~ sig_test$SEGMENT_GROUP)

## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/segment_groups_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/segment_groups_sig_assoc_POS_OR.csv", row.names = FALSE)

######
# Table 3 - Fisher's test domain groups
######

sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','DOMAIN_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")

# Loop through the domain groups
for (i in 1:length(domain_groups)) {
  
  dg = domain_groups[i]
  
  merged_1 <- merged %>% filter(!is.na(variant_type_2))
  
  # In domain group
  dg_terms <- merged_1 %>% filter(domain == dg)
  #pats = dg_terms$famID %>% unique()
  pats = dg_terms %>% filter(variant_type_2 == 'missense') %>% select(famID) %>% unique()
  pats <- as.vector(pats[,1])
  yes_hpo = merged_1 %>% filter(famID %in% pats)
  
  # Not in domain group
  no_hpo <- merged_1 %>% filter(variant_type_2 == 'missense')%>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','DOMAIN_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_cohort', 'absent_cohort', "group_freq", "cohort_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$DOMAIN_GROUP <- dg
  
  # Loop through hpos
  for(h in 1:nrow(hpo_1)){
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    # With var_x and hpo
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    # Without var_x and with hpo
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    # With var_x and without hpo
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    # Without var_x and hpo
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Raw counts and Fisher's test
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_cohort[h] = fish['group_absent','hpo_present']
    hpo_1$absent_cohort[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

sig_test <- sig_test %>% left_join(hpo_def)
sig_test <- sig_test %>% select(1, 14, 2:13)
sig_test$DOMAIN_GROUP <- case_when(TRUE~ sig_test$DOMAIN_GROUP)

## FDR

hpo_sig <- sig_test
fdr_adjust_or = FALSE # do not filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/domain_groups_sig_assoc_ALL_OR.csv", row.names = FALSE)

fdr_adjust_or = TRUE # filter or > 1
source("/Volumes/helbig_lab/projects/SCN2A/v13/primary_analyses/FDR.R")
write.csv(fdr_res, "frequency_analyses/domain_groups_sig_assoc_POS_OR.csv", row.names = FALSE)

