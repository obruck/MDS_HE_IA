# Created: OB, 01.12.2019

# This script runs a correlation analysis between TMA-level mutation label (predicted) and mutation variant allele frequency
## 0. Defines parameters --> lines 0-61
## 1. Recalls best-performing model to predict individual mutations --> lines 62-149
## 2. Joins mutation variant allele frequency (VAF) --> lines 150-197
## 3. Correlates predicted mutation probability with its actual VAF --> lines 198-258
## 4. Plot predicted mutation probability with its actual VAF for individual mutations --> lines 259-330

rm(list=ls())
# Load packages
library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
library(ggpubr)


##################################### Parameters ##############################################################################################


# Choose which parameter to run

# 1. Parameter "Y"
y_list = c("VGG16", "Xception")


# 2. Parameter "W" = what variables to run
## List of prediction labels is imported from a directory list
pred_filenames <- "/data/processed_files/Elastic_net_results"
file_list <- list.files(pred_filenames, full.names = FALSE, recursive = FALSE) 
### Remove unnecessary elements
w_list <- purrr::discard(file_list,.p = ~stringr::str_detect(.x,"November2019"))


# 3. Joining parameter z = what parameters to test for ROC
z_list <- purrr::discard(file_list,.p = ~stringr::str_detect(.x,"November2019")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"age_dg")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"Cmplx")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"cytopeniaipssr")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"IPSSR"))
z_list <- c("aza_log")


# 4. Parameter "regression type"
w_list_lin <- c("age_dg", "cytopeniaipssr","IPSSR")
w_list_cox <- c("aml", "overall_survival", "aza")
w_list_log <- purrr::discard(w_list,.p = ~stringr::str_detect(.x,"age_dg")) %>%
  purrr::discard(w_list,.p = ~stringr::str_detect(.x,"cytopeniaipssr")) %>%
  purrr::discard(w_list,.p = ~stringr::str_detect(.x,"IPSSR"))
w_list_log <- w_list_log[!w_list_log %in% w_list_lin]
w_list_log <- w_list_log[!w_list_log %in% w_list_cox]


# 5. Define alpha and lambda parameters
j_list = c("lasso_lambdamin", "lasso_lambda1se", "ridge_lambdamin", "ridge_lambda1se", "enet_lambdamin", "enet_lambda1se")


##################################### Best prediction models ##############################################################################################


# Define variables
w_list <- purrr::discard(file_list,.p = ~stringr::str_detect(.x,"November2019")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"aml_log")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"overall_survival_log")) %>%
  purrr::discard(file_list,.p = ~stringr::str_detect(.x,"aza"))
w_list <- w_list[!w_list %in% w_list_lin]
## Make lists for loop
auc_files = list()
auc_files_rownames = list()
auc_files_param_y = list()
auc_files_param_alpha = list()
auc_files_param_lambda = list()
## Generate AUC dataframe
for (y in y_list) {
  for (w in w_list) {
    for (j in j_list) {
      auc_files = rbind(auc_files, read_xlsx(paste0("/data/supervised_learning/", w, "/ROC/", y, "/ROC_", w, "_", y, "_", j, ".xlsx", sep="")))
      auc_files_param_y = c(auc_files_param_y, paste0(w))
      auc_files_param_alpha = c(auc_files_param_alpha, paste0(y))
      auc_files_param_lambda = c(auc_files_param_lambda, paste0(j))
      auc_files_rownames = c(auc_files_rownames, paste0("ROC_", w, "_", y, "_", j))
    }
  }
}


# Modify dataframe
auc_files <- auc_files %>%
  data.frame() %>%
  mutate(rowname = unlist(auc_files_rownames),
         auc_files_param_y = unlist(auc_files_param_y),
         auc_files_param_alpha = unlist(auc_files_param_alpha),
         auc_files_param_lambda = unlist(auc_files_param_lambda)) %>%
  column_to_rownames()
# Remove results of linear regressions
auc_files <- auc_files %>% filter(!auc_files_param_y %in% w_list_lin)
# Select based parameters based on test_wsi > train_wsi > test_tiles > train_tiles > enet>lasso>ridge
auc_files_best <- auc_files %>% 
  group_by(auc_files_param_y) %>% 
  arrange(-(test_wsi_auc), -(train_wsi_auc), -(test_tile_auc), -(train_tile_auc), auc_files_param_lambda) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(-(test_wsi_auc), -(train_wsi_auc), -(test_tile_auc), -(train_tile_auc), auc_files_param_lambda) %>%
  data.frame() %>%
  dplyr::select(5:7)


##################################### Mutation prediction models ##############################################################################################


# List of files to read
auc_files_best <- auc_files_best %>%
  mutate(input = paste0(auc_files_param_y, "/imagenet_", auc_files_param_alpha, "_", auc_files_param_lambda))
input <- as.list(auc_files_best$input)
## Remove non-mutation inputs
input <- purrr::discard(input,.p = ~stringr::str_detect(.x,"ch3abn")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"monosomy7")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"del7q")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"del5q")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"del20q")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"gender")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"abnchr")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"dg\\/")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"trisomy8")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"aml")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"Cmplx")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"overall_survival")) %>%
  purrr::discard(input,.p = ~stringr::str_detect(.x,"secondary_MDS"))


# Load prediction data for each mutation and each TMA image
# Dataframe structure
imagenet_cnn <- readRDS(paste0("/data/processed_files/Elastic_net_results/tet2/models/imagenet_Xception_lasso_lambda1se.rds"))
cnn_modellist <- as.data.frame(imagenet_cnn[4])[0,1:ncol(as.data.frame(imagenet_cnn[4]))]
## Load data
for (i in input) {
  imagenet_cnn <- readRDS(paste0("/data/processed_files/Elastic_net_results/", i, ".rds"))
  train_test_summary <- rbind(as.data.frame(imagenet_cnn[4]), as.data.frame(imagenet_cnn[5]))
  train_test_summary$input <- i
  train_test_summary$mut <- gsub("\\/[[:print:]]*", "", train_test_summary$input)
  cnn_modellist <- rbind(cnn_modellist, train_test_summary)
  imagenet_cnn <- NULL
}


##################################### Join mutation data ##############################################################################################


# Load mutation data
mut <- read_xlsx("/data/tables/MDS_mutations_table.xlsx", sheet = "Mutations_combined_vafs_unique")

# VAFs as numeric
mut[,which(colnames(mut) == "ASXL1"):ncol(mut)] <- sapply(mut[,which(colnames(mut) == "ASXL1"):ncol(mut)], function(x) as.numeric(as.character(x)))

# Select essential data
mut[is.na(mut)] = 0
mut <- mut %>%
  mutate(idh_mut = ifelse(IDH1 > IDH2, IDH1, IDH2)) %>%
  dplyr::select(patient_id, TET2, ASXL1, STAG2, SRSF2, TP53, RUNX1, SF3B1, ras_mut, idh_mut, Splicing, contains("mut_count"))

# Load metadata
metadata <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx") %>%
  dplyr::select(patient_id, Image_name)  %>% 
  rename(imagename = Image_name) %>%
  mutate(patient_id = as.numeric(patient_id)) %>%
  filter(patient_id > 20 & !is.na(patient_id))

# Merge metadata
mut <- mut %>% left_join(metadata)

# Merge prediction data
cnn_modellist1 <- mut %>% left_join(cnn_modellist)

# Rename mutations
cnn_modellist1 <- cnn_modellist1 %>%
  mutate(mut = ifelse(mut == "tet2", "TET2",
                      ifelse(mut == "asxl1", "ASXL1",
                             ifelse(mut == "stag2", "STAG2",
                                    ifelse(mut == "srsf2", "SRSF2",
                                           ifelse(mut == "tp53", "TP53",
                                                  ifelse(mut == "runx1", "RUNX1",
                                                         ifelse(mut == "sf3b1", "SF3B1",
                                                                ifelse(mut == "splicing_mut", "Splicing",
                                                                       ifelse(mut == "Celldifferentiation_mut", "Celldifferentiation_mut_count",
                                                                              ifelse(mut == "Cellcycle_mut", "Cellcycle_mut_count", 
                                                                                     ifelse(mut == "DNA_chromatin_mut", "DNA_chromatin_mut_count", 
                                                                                            ifelse(mut == "RAS_pathway_mut", "RAS_pathway_mut_count", mut)
                                                                                     ))))))))))))
# Filter samples with insufficient amount of cells for image analysis
cnn_modellist1 <- cnn_modellist1 %>% filter(!is.na(mut))


##################################### Correlate mutation vafs with prediction at TMA image level ##############################################################################################


# Correlate predictions with VAFs with all genes
mutations = unique(cnn_modellist1$mut)
cor_result = data.frame(matrix(ncol = 3, nrow = 0))
cor_result = list()
for (i in mutations) {
  tmp <- cnn_modellist1 %>% dplyr::filter(mut == i)
  tmp <- tmp[tmp[,i] > 0,]
  tmp <- as.data.frame(tmp)
  cor_summary <- cor.test(tmp[,i], tmp$image_pred_mean, method = "spearman", exact = FALSE)
  cor_result1 <- list(c(i, cor_summary$estimate, cor_summary$p.value))
  cor_result <- c(cor_result, cor_result1)
}


## Convert list of lists by dataframe
cor_result <- do.call(rbind.data.frame, cor_result)
colnames(cor_result) <- c("mutation", "cor", "pvalue")


## Rename mutations by adding "_mut" suffix
cor_result <- cor_result %>% mutate(
  mutation = gsub("2$", "2 VAF", mutation),
  mutation = gsub("1$", "1 VAF", mutation),
  mutation = gsub("3$", "3 VAF", mutation),
  mutation = gsub("RAS_pathway_mut_count", "RAS pathway mutation count", mutation),
  mutation = gsub("DNA_chromatin_mut_count", "Chromatin pathway mutation count", mutation),
  mutation = gsub("Cellcycle_mut_count", "Cell cycle pathway mutation count", mutation),
  mutation = gsub("Celldifferentiation_mut_count", "Cell differentiation pathway mutation count", mutation),
  mutation = gsub("Splicing", "Spliceosome mutation VAF", mutation),
  mutation = gsub("ras_mut", "KRAS/NRAS VAF", mutation),
  mutation = gsub("idh_mut", "IDH1/IDH2 VAF", mutation)
) %>% mutate(
  cor = as.numeric(as.character(cor)),
  pvalue = as.numeric(as.character(pvalue))
)


# Plot
png("/data/results/Correlation_scatter_plot.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(cor_result, aes(x=cor, y=-log10(pvalue), fill = mutation, size = cor)) +
  geom_point(colour="black", pch=21) +
  geom_hline(aes(lty = "p<0.05", yintercept=10^0.05)) +
  labs(
    size = "Cor",
    fill = "Mutation",
    linetype = "P",
    x="Spearman Correlation Coefficient",
    y="-Log10 P-value") +
  theme(axis.title.y = element_text(size=12,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        legend.position="bottom") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#fdbf6f", "#b2df8a", "black", "#cab2d6", "#e7298a")) +
  guides(size = guide_legend(ncol = 1),
         linetype = FALSE,
         fill = guide_legend(ncol = 3, override.aes = list(size = 5))) +
  scale_linetype_manual(values=2)
dev.off()


##################################### Export images ##############################################################################################


# Produce correlation plots
mutations = unique(cnn_modellist1$mut)
mutations_cor <- purrr::discard(mutations,.p = ~stringr::str_detect(.x,"mut_count"))
mutations_cor <- mutations
for (i in mutations_cor) {
  tmp <- cnn_modellist1 %>% dplyr::filter(mut == i)
  tmp <- tmp[tmp[,i] > 0,]
  tmp <- as.data.frame(tmp)
  
  
  ## Rename mutations by adding "_mut" suffix
  j <- gsub("2$", "2 mutation", i)
  j <- gsub("1$", "1 mutation", j)
  j <- gsub("3$", "3 mutation", j)
  j <- gsub("RAS_pathway_mut_count", "mutations in RAS pathway", j)
  j <- gsub("DNA_chromatin_mut_count", "mutations in DNA chromatin regulation pathway", j)
  j <- gsub("Cellcycle_mut_count", "mutations in cell cycle regulation pathway", j)
  j <- gsub("Celldifferentiation_mut_count", "mutations in cell differentiation regulation pathway", j)
  j <- gsub("Splicing", "splicing mutation", j)
  j <- gsub("ras_mut", "RAS pathway mutation", j)
  j <- gsub("idh_mut", "IDH pathway mutation", j)
  
  
  g <- ggplot(tmp, aes_string(x="image_pred_mean", y=i)) +
    geom_point(size=4) +
    geom_smooth(method = "lm", color="black") +
    labs(
      x=paste0("Proportion of ", j),
      y=paste0("Prediction probability of ", j)) +
    theme(axis.title.y = element_text(size=12,face="bold"),
          axis.title.x = element_text(size=12,face="bold")) +
    scale_color_manual("black") +
    stat_cor(method = "spearman", size = 6)
  ggsave(g, file=paste0("/data/results/Correlation_", i, ".png"), width = 6, height = 6, limitsize = FALSE, units = 'in', dpi = 300)
}

# Produce correlation plots
mutations = unique(cnn_modellist1$mut)
mutations_cor <- purrr::keep(mutations,.p = ~stringr::str_detect(.x,"mut_count"))
for (i in mutations_cor) {
  tmp <- cnn_modellist1 %>% dplyr::filter(mut == i)
  tmp <- tmp[tmp[,i] > 0,]
  tmp <- as.data.frame(tmp)
  
  
  ## Rename mutations by adding "_mut" suffix
  j <- gsub("2$", "2 mutation", i)
  j <- gsub("1$", "1 mutation", j)
  j <- gsub("3$", "3 mutation", j)
  j <- gsub("RAS_pathway_mut_count", "mutations in RAS pathway", j)
  j <- gsub("DNA_chromatin_mut_count", "mutations in DNA chromatin regulation pathway", j)
  j <- gsub("Cellcycle_mut_count", "mutations in cell cycle regulation pathway", j)
  j <- gsub("Celldifferentiation_mut_count", "mutations in cell differentiation regulation pathway", j)
  j <- gsub("Splicing", "splicing mutation", j)
  j <- gsub("ras_mut", "RAS pathway mutation", j)
  j <- gsub("idh_mut", "IDH pathway mutation", j)
  
  
  g <- ggplot(tmp, aes_string(x="image_pred_mean", y=i)) +
    geom_point(size=4) +
    labs(
      x=paste0("Proportion of ", j),
      y=paste0("Prediction probability of ", j)) +
    theme(axis.title.y = element_text(size=12,face="bold"),
          axis.title.x = element_text(size=12,face="bold")) +
    scale_color_manual("black") +
    stat_cor(method = "spearman", size = 6)
  ggsave(g, file=paste0("/data/results/Correlation_", i, "_2.png"), width = 6, height = 6, limitsize = FALSE, units = 'in', dpi = 300)
}

