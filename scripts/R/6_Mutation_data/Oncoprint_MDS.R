# Created: 18.10.2019, Oscar Brück

# Oncoprint

rm(list=ls())

# Load packages
library(tidyverse)
library(lubridate)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


######################## Load data ############################################################################################################################


# Mutations
mut <- read_xlsx("/data/tables/MDS_mutations_table.xlsx", sheet = "Mutations_combined")
## Only mutation data from HUS
mut <- mut %>% dplyr::select(c(patient_id, Batch, 7:(ncol(mut)-10))) %>% dplyr::select(-c(contains("_KI"), EP300_HKI))

# Mutation groups
mut_groups <- read_xlsx("/data/tables/Gene_groups.xlsx", sheet = "Vertical")

# Diagnoses
mut_dg <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx") %>%
  mutate( patient_id = as.numeric(patient_id) ) %>%
  arrange(patient_id, Saapunut) %>%
  filter(!is.na(patient_id)) %>%
  filter(diagnosis_abbr_sample_HEproject == "MDS" | diagnosis_abbr_sample_HEproject == "PreMDS") %>%
  mutate(diagnosis_pt = ifelse(diagnosis_abbr_sample_HEproject == "MDS" | diagnosis_abbr_sample_HEproject == "PreMDS", "MDS", NA)) %>%
  filter(!is.na(diagnosis_pt)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  dplyr::select(patient_id, diagnosis_pt)


#########################################################################################################################################################


# Replace
## Move first two columns
mut1 <- mut[,1:2]
mut[,1:2] <- NULL

## Mutations to MUT
mut[!is.na(mut) & mut != "UNKNOWN"] <- "MUT"
## NA with ""
mut[is.na(mut)] <- ""
## Return first two columns
mut <- cbind(mut1, mut)


# 10 patients have been sequenced twice --> summarise
mut <- mut %>% dplyr::select(-Batch) %>% group_by(patient_id) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ''))
mut <- sapply(mut, function(x) gsub("UNKNOWNMUT", "MUT", x))
mut <- data.frame(mut)
mut <- sapply(mut, function(x) gsub("MUTUNKNOWN", "MUT", x))
mut <- data.frame(mut)
mut$patient_id <- as.numeric(as.character(mut$patient_id))

# Add diagnosis information
mut <- mut %>% left_join(mut_dg) %>% filter(!is.na(diagnosis_pt))


# As dataframe and matrix
mut1 <- NULL
mut1 <- mut %>% dplyr::select(-patient_id, -diagnosis_pt)

# Column order
colnames(mut1) <- gsub("_KI", "", colnames(mut1))
colnames(mut1) <- gsub("_HKI", "", colnames(mut1))
mut1_summary_pergene <- sapply(mut1, function(x) sum(x == "MUT"))
mut1_summary_pergene <- mut1_summary_pergene[mut1_summary_pergene>2]


# Filter genes with frequency < 3
mut1 <- mut1[,colnames(mut1) %in% names(mut1_summary_pergene)]


# Row order
mut1_summary_perpt <- sapply(as.data.frame(t(mut1)), function(x) sum(x == "MUT"))


# Gene groups
mut_groups <- mut_groups %>% filter(Gene %in% colnames(mut1)) %>% dplyr::select(-Leukemia)
tmp <- data.frame(mut1_summary_pergene) %>% rownames_to_column() %>% rename(Gene = rowname)
mut_groups <- mut_groups %>% dplyr::select(Gene, Genegroup_tmp2_2) %>% left_join(tmp)

## Merge
mut1 <- mut1 %>% t() %>% data.frame() %>% rownames_to_column() %>% rename(Gene = rowname)
mut1 <- mut1 %>% left_join(mut_groups)
mut1 <- mut1 %>% t() %>% data.frame() 
colnames(mut1) <- unlist(mut1[1,])
mut1 <- mut1[2:nrow(mut1),]
rownames(mut1) <- NULL

# Bind mut1_summary_pergene with mut1
mut2 <- as.matrix(t(mut1))


####################################### Oncoprint #############################################################################################


# Oncoprint grid
# #C77CFF, , #7CAE00
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = "black"))  # col = NA
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F8766D", col = "black"))
  },
  UNKNOWN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#00BFC4", col = "black"))
  }
)

col = c("MUT" = "#F8766D", "background" = "#e0e0e0")



# Annotation bars

# Top annotation
top_annotation <- HeatmapAnnotation(
  Freq = anno_barplot(mut1_summary_perpt, axis = TRUE, gp=gpar(border=NA, fill="#F8766D", col = "black"),
                      bar_width = 1))
  
# Left annotation
left_annotation <- rowAnnotation("%" = anno_text(paste(round(100 * as.numeric(mut2[,ncol(mut2)]) / ncol(mut2), 0), "%", sep="")))


# Oncoprint plot
ht_list1 <- oncoPrint(mut2[, 1:(ncol(mut2)-2)],
          alter_fun = alter_fun,
          col = col,
          show_pct = FALSE,
          left_annotation = left_annotation,
          row_order = order(-as.numeric(mut2[,ncol(mut2)])),
          remove_empty_columns = FALSE,
          remove_empty_rows = TRUE,
          top_annotation = top_annotation,
          column_title = "Mutations in Driver Genes of MDS Patients",
          heatmap_legend_param = list(title = "Legend", at = c("MUT"), labels = c("Mutation"))
)
ht_list1



ht_list2 <- Heatmap(mut2[,(ncol(mut2)-1)],
                    name = "Gene Group",
                    rect_gp = gpar(col = "black"),
                    show_column_dend = FALSE,
                    col = c("DNAchromatin" = "#e41a1c",
                            "Cellcycle" = "#377eb8",
                            "Celldifferentiation" = "#4daf4a",
                            "Splicing" = "#a65628",
                            "Cellcycle&Celldifferentiation" = "#984ea3",
                            "DNAchromatin&Cellcycle" = "#ff7f00",
                            "DNAchromatin&Cellcycle&Celldifferentiation" = "#ffff33"),
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    width = unit(0.5, "cm")
)


png("/data/results/Oncoprint_MDS.png", width = 11.5, height = 6, units = 'in', res = 300)
draw(ht_list1+ht_list2, newpage = TRUE, row_sub_title_side = "left", merge_legend = TRUE)
dev.off()


