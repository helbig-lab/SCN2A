library(tidyverse)
library(Hmisc)
library(logisticPCA)
library(reshape2)

####################################
# 1. PCA
####################################

################## 
# Data input
################## 

scn2a <- read.csv(input.yaml$variant_file) %>%
  filter(variant_type_1 != "exclude") %>%
  filter(grepl("HP:[0-9]", HPO)) %>%
  unique()

func_variants <- read.csv(input.yaml$functional_variant, stringsAsFactors = F) %>%
  rename(effect = Overall.Effect) %>%
  select(variant, effect)

scn2a_pos <- read.csv(input.yaml$pos_prop)
scn2a_neg <- read.csv(input.yaml$neg_prop)

scn2a_prop <- scn2a_pos %>% rbind(scn2a_neg)
rm(scn2a_pos, scn2a_neg)

# Hpo dictionary
pos_hpo <- read.csv(input.yaml$hpo_ancestor) %>%
  select(term, def) %>%
  rename(HPO = term) %>%
  mutate(HPO = trimws(HPO))

################## 
# Data cleaning
################## 

variants_all <- scn2a %>% 
  select(variant, variant_type_2) %>% 
  unique() %>% 
  filter(variant_type_2 != "") %>% 
  left_join(func_variants) %>% 
  mutate(functional = case_when(variant_type_2 == "missense" & effect == "GOF" ~ "Missense GOF",
                                variant_type_2 == "missense" & effect == "LOF" ~ "Missense LOF",
                                variant_type_2 == "missense" & effect %nin% c("GOF", "LOF") ~ "Missense other",
                                variant_type_2 == "PTV" ~ "PTV")) %>% 
  mutate(functional_1 = case_when(functional == "Missense GOF" ~ "Missense GOF",
                                  functional %in% c("Missense LOF", "PTV") ~ "Missense LOF + PTV",
                                  functional == "Missense other" ~ "Missense other"))

scn2a_prop <- scn2a_prop %>% 
  left_join(scn2a %>% select(famID, variant) %>% unique()) %>% 
  left_join(variants_all %>% select(variant, variant_type_2, functional_1) %>% unique())


# Add negatives to hpo dictionary
neg_hpo <- pos_hpo %>% 
  mutate(HPO = gsub("H","N",HPO)) %>% 
  mutate(def = paste("No ",def, sep = "")) %>% 
  unique()

hpo_all <- pos_hpo %>% rbind(neg_hpo)
rm(neg_hpo)

scn2a <- scn2a %>% filter(!is.na(variant_type_2))

# NROW(unique(scn2a$variant[scn2a$function.=='GOF']))

################## 
# Load pca
################## 

# Output from 'run_variant_ks.R'
load("functional_variant_logpca_model_k3.RData")

################## 
# Analysis
################## 

pats = as.vector(scn2a$famID %>% unique())
hpo = as.vector(scn2a_prop$HPO %>% unique())

# PC values
pca_x <- as.data.frame(logpca_model$PCs); pca_x$famID = pats
colnames(pca_x) = c("PC1", "PC2", "PC3", "famID")
pca_x <- pca_x %>% select(famID, 1:(ncol(pca_x)-1)) %>% left_join(scn2a_prop %>% select(-HPO) %>% unique())

# Variance explained
logpca_model$prop_deviance_expl
vars_transformed <- apply(pca_x %>% select("PC1", "PC2", "PC3"), 2, var)
vars_transformed/sum(vars_transformed)*logpca_model$prop_deviance_expl

# Factor loadings
loadings <- as.data.frame(logpca_model$U); loadings$HPO = hpo
colnames(loadings) = c("PC1", "PC2", "PC3", "HPO")
loadings <- loadings %>% left_join(hpo_all) %>% select(HPO, def, 1:(ncol(loadings)-1))

###

####################################
# 2. ROC
####################################

#https://hopstat.wordpress.com/2014/12/19/a-small-introduction-to-the-rocr-package/

library(tidyverse)
library(ROCR)
library(corrplot)
library(RColorBrewer)
library(Hmisc)


#########
# Data input
#########

# Requires PC values

roc_x <- pca_x
rm(list=setdiff(ls(), "roc_x"))

# Get "optimal" cut point - minimum distance from (0,1)
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

#########
# ROC for Figure
#########

# PC = 2
roc <- list()
auc_all <- c()
f1_all <- c()
pc_aggregate <- list() #MCK edit - to be used later for multiple PR curves
prec_recall <- list() # JX
for (pcx in 1:3) {
  
  PC = pcx
  
  g1 = "Missense LOF + PTV"
  g2 = "Missense GOF"
  
  # Decide what group to use 0 and 1
  df <- roc_x %>% 
    select(famID, paste0("PC", PC), functional_1) %>% 
    filter(functional_1 %in% c(g1, g2)) 
  
  if (pcx == 2) {
    df_PC2 = df
  }
  
  m1 = median(df %>% filter(functional_1 == g1) %>% select(paste0("PC", PC)) %>% pull())
  m2 = median(df %>% filter(functional_1 == g2) %>% select(paste0("PC", PC)) %>% pull())
  
  if (m1 < m2) {
    df <- df %>% mutate(pred = case_when(functional_1 == g1 ~ 0, functional_1 == g2 ~ 1))
  } else {
    df <- df %>% mutate(pred = case_when(functional_1 == g1 ~ 1, functional_1 == g2 ~ 0))
  }
  
  pc_aggregate[[pcx]] <- df #MCK Edit - to be used for PR later
  
  pred_std <- prediction(df %>% select(paste0("PC", PC)), df$pred)
  perf_std <- performance(pred_std,"tpr","fpr")
  
  #MCK Addition - F1 scores calculated and added to cutoff_std
  f1_std <- performance(pred_std, measure="f")
  
  cutoffs_std <- data.frame(cutoff=perf_std@alpha.values[[1]], 
                            fpr=perf_std@x.values[[1]], 
                            tpr=perf_std@y.values[[1]],
                            f1=f1_std@y.values[[1]])
  
  # Calculate PPV for PC2
  # if (pcx == 2) {
  
  ppv_std <- cutoffs_std
  ppv_std$ppv_gof = NA
  ppv_std$ppv_lof = NA
  for (c in 1:(nrow(ppv_std))) {
    
    cut <- ppv_std$cutoff[c]
    
    # cut_l <- df %>% filter(PC2 <= cut)
    # cut_r <- df %>% filter(PC2 > cut)
    
    cut_l <- df[df[,2] <= cut,]
    cut_r <- df[df[,2] > cut,]
    
    #At loading cutoff, percent GOF and percent LOF
    if (nrow(cut_r) > 0 & nrow(cut_l) > 0) {
      ppv_std$ppv_gof[c] <- sum(cut_r == 1)/(sum(cut_r == 1) + sum(cut_r == 0))
      ppv_std$ppv_lof[c] <- sum(cut_l == 0)/(sum(cut_l == 1) + sum(cut_l == 0))
    } else if (nrow(cut_l) == 0) {
      ppv_std$ppv_gof[c] <- sum(cut_r == 1)/(sum(cut_r == 1) + sum(cut_r == 0))
      ppv_std$ppv_lof[c] <- NA
    } else if (nrow(cut_r) == 0) {
      ppv_std$ppv_lof[c] <- sum(cut_l == 0)/(sum(cut_l == 1) + sum(cut_l == 0))
      ppv_std$ppv_gof[c] <- NA
    }
    
  }
  ppv_std <- ppv_std %>% na.exclude()
  
  prec_recall[[pcx]] <- ppv_std # JX
  
  ppv_df <- ppv_std %>%
    select(cutoff, ppv_lof) %>%
    rename(ppv = ppv_lof) %>%
    mutate(pred = "Missense LOF + PTV") %>%
    rbind(ppv_std %>%
            select(cutoff, ppv_gof) %>%
            rename(ppv = ppv_gof) %>%
            mutate(pred = "Missense GOF"))
  # }
  
  cutoffs_std <- cutoffs_std #%>% select(-cutoff) MCK Edit
  
  
  
  roc[[pcx]] <- cutoffs_std %>% select(-cutoff)
  
  #MCK Edit
  opt_std = opt.cut(perf_std, pred_std);print(opt_std); sen=opt_std[1]; spec=opt_std[2]
  
  print("F1")
  opt_f1 = cutoffs_std$f1[cutoffs_std$cutoff == opt_std[3]] ;print(opt_f1)
  
  # Save cutoff for PC2
  if (pcx == 2) {
    opt_cut_PC2 = opt_std[3]
  }
  
  # Get the area under the curve (AUC)
  auc.perf_std = performance(pred_std, measure = "auc"); auc = auc.perf_std@y.values[[1]]
  x=auc
  print(paste0("AUC ",x))
  auc_all <- append(auc_all, x)
  
}

# ROC plot (PC 1-3)

df.long <- melt(roc,id.vars ="fpr") %>% rename(PC = L1) %>% mutate(PC = paste0("PC",PC))
df.long <- df.long %>% mutate(AUC = case_when(PC == "PC1" ~ paste0("PC1 (AUC = ", round(auc_all[1],2), ")"),
                                              PC == "PC2" ~ paste0("PC2 (AUC = ", round(auc_all[2],2), ")"),
                                              PC == "PC3" ~ paste0("PC3 (AUC = ", round(auc_all[3],2), ")")))


## CHANGE SCALE COLOR MANUAL

png(filename=paste0("variant_ROC.png"),
    units = "in", width = 10, height = 10, res = 1000)

(pcp <- ggplot(df.long %>% filter(PC %in% c("PC1", "PC2", "PC3")),aes(fpr,value,color=AUC)) + 
    geom_step(size = 0.75) +
    scale_color_manual(values = c("PC1 (AUC = 0.57)" = "steelblue",
                                  "PC2 (AUC = 0.84)" = "orange",
                                  "PC3 (AUC = 0.63)" = "slategray")) +
    labs(x="False positive rate",y="True positive rate") +
    geom_abline(intercept=0, slope=1, linetype="dashed", col="gray") +
    theme_classic() +
    theme(legend.position = c(.8,.15), legend.title = element_blank(), legend.text=element_text(size=16),
          axis.text = element_text(size = 12), axis.title = element_text(size = 16), legend.key.height = unit(2,"line")))


dev.off()


###

#MCK!!!
# PR Curve Plot - Take each PC dataframe from above and reshape as precision/recall datasets for each PC
pr_list = list()
for (i in 1:NROW(pc_aggregate)) {
  
  pred_std <- prediction(as.data.frame(pc_aggregate[i]) %>% select(paste0("PC", i)), as.data.frame(pc_aggregate[i]) %>% select(pred)) # JX changed
  perf_pr <- performance(pred_std, "prec", "rec")
  
  perf_pr
  plot(perf_pr)
  
  pr_df <- data.frame(rec=perf_pr@x.values,
                      prec=perf_pr@y.values, 
                      cutoff=perf_pr@alpha.values, PC=paste0("PC", i))
  colnames(pr_df) <-  c("recall", "precision", "cutoff", "PC")
  
  pr_list[[i]] <- pr_df
  
}
pr_dataframe = do.call(rbind, pr_list)

ggplot(pr_dataframe, aes(x=recall, y=precision, col=PC)) +
  geom_step(size = 0.75) +
  #geom_abline(intercept = 0) +
  labs(x = "Recall", y = "Precision") +
  scale_color_manual(values = c("PC1" = "steelblue",
                                "PC2" = "orange",
                                "PC3" = "slategray")) +
  scale_x_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) +
  #coord_cartesian(xlim = c(0, 1), ylim = (0, 1)) +
  #geom_text(aes(x = .5, y = .5, label = PC, color = PC), check_overlap = T) + 
  theme_classic() +
  theme(legend.position = c(.8,.15), legend.title = element_blank(), legend.text=element_text(size=16),
        axis.text = element_text(size = 12), axis.title = element_text(size = 16), legend.key.height = unit(2,"line"))
#ggsave(filename=paste0("/Volumes/helbig_lab/projects/SCN2A/v16/figures/Fig_S9/variant_PR.png"),
#            units = "in", width = 10, height = 10, dpi = 1000)

dev.off()





# PPV plot


variants <- roc_x %>% select(famID, PC2, variant, variant_type_2, functional_1) %>% 
  left_join(scn2a %>% select(famID, broad_phx) %>% unique())

variants$cutoff = NA
for (i in 1:NROW(variants)) {
  variants$cutoff[i] <- ppv_df[which(abs(ppv_df$cutoff-variants$PC2[i]) == min(abs(ppv_df$cutoff-variants$PC2[i]))),]$cutoff %>% unique()
}
ppv_variants <- variants %>% left_join(ppv_df)

## Table
ppv_variants_1 <- variants %>% 
  left_join(ppv_std) %>% 
  select(famID, variant, functional_1, PC2, ppv_gof, ppv_lof)

count <- ppv_variants_1 %>% 
  dplyr::group_by(variant) %>% 
  count()

ppv_variants_1 <- ppv_variants_1 %>% 
  left_join(count) %>% 
  rename(variant_n = n) %>% 
  select(1,2,7,3,4,5,6)

ppv_variants_1[ppv_variants_1$variant == "E1211K",]$functional_1 = "Mixed"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense GOF",]$functional_1 = "GoF"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense LOF + PTV",]$functional_1 = "LoF"
ppv_variants_1[ppv_variants_1$functional_1 == "Missense other",]$functional_1 = "none"

write.csv(ppv_variants_1, "fam_PC2_PPV.csv", row.names = F)

##


ppv_variants <- ppv_variants %>% mutate(var_lab = case_when(famID %in% c("fam34", "fam94", "fam349", "fam92", 
                                                                         "fam250", "fam83", "fam196", "fam235") 
                                                            & pred == "Missense GOF" ~ "yes",
                                                            famID %in% c("fam302", "fam112", "fam108", "fam144", 
                                                                         "fam52", "fam362","fam355", "fam163", "fam226") 
                                                            & pred == "Missense LOF + PTV" ~ "yes", 
                                                            TRUE ~ "no"))

png(filename=paste0("PPV_PC2_1.png"),
    units = "in", width = 10, height = 10, res = 1000)

ggplot(ppv_variants, aes(x = cutoff, y = ppv, color = pred)) +
  geom_line(size = 0.75) +
  theme_bw() +
  labs(x = "PC 2", y = "Positive predictive value") +
  scale_color_manual(values = c("Missense GOF" = "red",
                                "Missense LOF + PTV" = "blue",
                                "Missense other" = "grey60")) +
  ggrepel::geom_label_repel(aes(label = ifelse(var_lab == "yes", as.character(variant),'')), 
                            hjust=0, vjust=0,colour = "black", size = 7, min.segment.length = 0.5) +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16))

dev.off()

# Without labels
png(filename=paste0("PPV_PC2_nl_1.png"),
    units = "in", width = 10, height = 10, res = 1000)

ggplot(ppv_variants, aes(x = cutoff, y = ppv, color = pred)) +
  geom_line(size = 0.75) +
  theme_bw() +
  labs(x = "PC 2", y = "Positive predictive value") +
  scale_color_manual(values = c("Missense GOF" = "red",
                                "Missense LOF + PTV" = "blue",
                                "Missense other" = "grey60")) +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16))


dev.off()


### Step plot

png(filename=paste0("PPV_PC2_step_1.png"),
    units = "in", width = 10, height = 10, res = 1000)

ggplot(ppv_variants, aes(x = cutoff, y = ppv, color = pred)) +
  geom_step(size = 0.75) +
  theme_bw() +
  labs(x = "PC 2", y = "Positive predictive value") +
  scale_color_manual(values = c("Missense GOF" = "red",
                                "Missense LOF + PTV" = "blue",
                                "Missense other" = "grey60")) +
  ggrepel::geom_label_repel(aes(label = ifelse(var_lab == "yes", as.character(variant),'')), 
                            hjust=0, vjust=0,colour = "black", size = 7, min.segment.length = 0.5) +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16))

dev.off()

# Without labels
png(filename=paste0("PPV_PC2_nl_step_1.png"),
    units = "in", width = 10, height = 10, res = 1000)

ggplot(ppv_variants, aes(x = cutoff, y = ppv, color = pred)) +
  geom_step(size = 0.75) +
  theme_bw() +
  labs(x = "PC 2", y = "Positive predictive value") +
  scale_color_manual(values = c("Missense GOF" = "red",
                                "Missense LOF + PTV" = "blue",
                                "Missense other" = "grey60")) +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16))

dev.off()



###
###



# Standard PCA plot
(p <- ggplot(roc_x, aes(x = PC2, y = PC1)) +
    geom_point(aes(colour = functional_1)) +
    theme_minimal() +
    scale_color_manual(values = c("Missense GOF" = "red",
                                  "Missense LOF + PTV" = "blue",
                                  "Missense other" = "grey80")) +
    theme(legend.position = "none"))



# Only PC 2
(p_o2 <- ggplot(roc_x) +
    theme_classic() +
    aes(x = roc_x %>% select("PC2") %>% pull(), y = paste0("PC2")) +
    geom_point(alpha = 0.8, position = position_jitter(width = 0.1), aes(color = functional_1, size = 2)) +
    # geom_boxplot(alpha = 0, size = 0.75, width = 0.25) +
    theme(legend.position = "none", legend.title = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_line(colour = "gray", size = 0)) +
    scale_color_manual(values = c("Missense GOF" = "red",
                                  "Missense LOF + PTV" = "blue",
                                  "Missense other" = "grey80")) +
    labs(y="", x = "") +
    geom_vline(xintercept=opt_cut_PC2, linetype="dashed", col="gray30"))

ggsave(plot = p_o2, filename = "PCA_PC2.png",dpi = 1000, width = 8, height = 5)


roc_x_1 <- roc_x %>% mutate(functional_1 = case_when(functional_1 == "Missense GOF" ~ "GoF",
                                                     functional_1 == "Missense LOF + PTV" ~ "LoF",
                                                     functional_1 == "Missense other" ~ "Other"))


# Density
(d <- ggplot(roc_x_1, aes(x = PC2, colour = functional_1)) +
    geom_density(position="identity", fill = NA, size = 1) +
    theme_classic() +
    scale_colour_manual(values = c("GoF" = "red",
                                   "LoF" = "blue",
                                   "Other" = "grey60")) +
    labs(y = "Density") +
    geom_vline(xintercept=opt_cut_PC2, linetype="dashed", col="gray30") +
    theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 16)))

ggsave(plot = d, filename = "Density_PC2_nf.png",dpi = 1000, width = 8, height = 5)

