# Created: 31.03.2019, OB
# Modified: 31.03.2019, OB

# MDS HE unsupervised learning - cluster statistics

rm(list=ls())

# Load packages
library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
library(ggimage)
library(ggplot2)
library(ggpubr)


############### Parameters ###############################################################################


# Define which CNN network dataset to use
y_list <- "Xception"

# Define which diseases to include
dg_list1 <- c("MDS", "MDS_MDSMPN", "MDS_MDSMPN_CO", "MDS_MDSMPN_AA_CO", "MDS_CO")

# Define which disease stage to analyze
dg_phase1 <- c("dg", "dg_fu")

# Define WSI or tile
image_size <- c("wsi")

## Clustering
clustering <- c("kcluster_umap", "phenograph_cluster_umap")


########################## Looping to calculate intercluster statistics ##################################################################


for (dg_list in dg_list1) {
  for (y in y_list) {
    for (dg_phase in dg_phase1) {
      print(paste0("Exporting statistics for ", dg_list, " ", dg_phase, " data, generated with the ", y, " framework"))
      # }}}

      
      ########################## Export and import data for plots ##################################################################
      
      
      # Change working directory
      setwd(paste0("/data/results/Unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"))
      print(getwd())
      
      
      # Export and import data for plots
      vgg16_feature_list_cluster <- readRDS(paste0("Unsupervised_learning_", dg_list, ".rds"))
      
      
      ########### Modify variables ##############################################################################
      
      
      print("Modify data")
      
      # Select labels for plots
      # Modify data
      vgg16_feature_list_cluster <- vgg16_feature_list_cluster %>% rename(
        Diagnosis = diagnosis_abbr_sample_HEproject,
        Stroma_proportion = stroma_proportion,
        Lipid_droplet_proportion = lipid_droplet_proportion_of_mask,
        RBC_proportion = rbc_proportion,
        WBC_proportion = tumor_proportion,
        IPSS_score = ipss_score,
        IPSS_class = ipss_class,
        IPSSR_score = ipssr_score,
        IPSSR_class = ipssr_class
      ) %>% mutate(
        Dysplasia = ifelse(Dysplasia_HUSLAB==1, "Dysplasia", "No dysplasia"),
        cmplx_ipss = ifelse(is.na(cmplx_ipss), "NA",
                            ifelse(cmplx_ipss==1, "Complex", "Non-complex")),
        Complex_karyotype = factor(cmplx_ipss, levels=c("Complex", "Non-complex", "NA")),
        ch_five_q_del = ifelse(is.na(ch_five_q_del), "NA",
                               ifelse(ch_five_q_del==1, "5qdel", "Non-5qdel")),
        Del5q = factor(ch_five_q_del, levels=c("5qdel", "Non-5qdel", "NA")),
        diagnose_kval = ifelse(is.na(diagnose_kval), "NA",
                               ifelse(diagnose_kval==1, "Abnormal karyotype", "Normal karyotype")),
        Karyotype = factor(diagnose_kval, levels=c("Abnormal karyotype", "Normal karyotype", "NA")),
        Disease_status = factor(vgg16_feature_list_cluster$Disease_status, levels=c("Dg", "FU"))
      ) %>% mutate(
        ### Disease status
        icd_o_pt_HEproject = ifelse(icd_o_pt_HEproject == "MDS/MPN + ring sideroblasts", "MDS/MPN-U",
                                    ifelse(icd_o_pt_HEproject == "9975 Myelodysplastic/Myeloproliferative neoplasm, unclassifiable", "MDS/MPN-U", 
                                           ifelse(icd_o_pt_HEproject == "Refractory cytopenia with singlelineage dysplasia", "MDS-SLD", 
                                                  ifelse(icd_o_pt_HEproject == "9985 Refractory cytopenia with multilineage dysplasia", "MDS-MLD", 
                                                         ifelse(icd_o_pt_HEproject == "9982 Refractory anemia with sideroblasts", "MDS-RS", 
                                                                ifelse(icd_o_pt_HEproject == "9945 Chronic myelomonocytic leukemia, NOS", "CMML", 
                                                                       ifelse(icd_o_pt_HEproject == "9986 Myelodysplastic syndr. with 5q deletion syndrome", "MDS-Del(5q)", 
                                                                              ifelse(icd_o_pt_HEproject == "MDS, unclassifiable (MDS-U)", "MDS-U", 
                                                                                     ifelse(Diagnosis == "MDS" & is.na(icd_o_pt_HEproject), "MDS-U", 
                                                                                            ifelse(Diagnosis == "MDS/MPN" & is.na(icd_o_pt_HEproject), "MDS/MPN-U",
                                                                                                   ifelse(Diagnosis == "MDS" & icd_o_pt_HEproject == "MDS/MPN-U", "MDS-U", 
                                                                                                          icd_o_pt_HEproject))))))))))),
        Diagnosis = ifelse(Diagnosis == "PreMDS", "MDS", Diagnosis),
        Diagnosis = as.factor(Diagnosis),
        IPSS_class = factor(IPSS_class, levels = c("Low", "Int-1", "Int-2", "High", NA)),
        IPSSR_class = factor(IPSSR_class, levels = c("Very low", "Low", "Intermediate", "High", "Very high", NA))
      ) %>%
        rename(WHO = icd_o_pt_HEproject)
      
      
      ############### Variable parameters #######################################################################
      
      
      print("Define parameters")
      
      # All variables
      mutations <- vgg16_feature_list_cluster %>% dplyr::select(contains("_kval"), -"diagnose_kval") %>% colnames()
      
      variable_names <- c("Diagnosis", "Dysplasia", "IPSS_class", "IPSSR_class", "Complex_karyotype", "Del5q", "Karyotype", "Disease_status", "WHO",
                          "Stroma_proportion", "Lipid_droplet_proportion", "RBC_proportion", "WBC_proportion", "IPSS_score", "IPSSR_score", "secondary_MDS",
                          "Ymonosomy", "Del11q", "Del12p", "Del20q", "Del7q", "Monosomy7", "Trisomy8", "Trisomy19", "i17q", "Ch3abn",
                          "cytogenetics_ipss_class", "cytogenetics_ipssr_class", "cytogenetics_ipss_class_score", "cytogenetics_ipssr_class_score",
                          "b-hb", "b-leuk", "b-neut", "b_trom", "bm_blast_p", "cytopenia_score_ipss", "cytopenia_score_ipssr",
                          "Qupath_Nucleus_Area", "Qupath_Nucleus_Perimeter", "Qupath_Nucleus_Circularity", "Qupath_Nucleus_Eccentricity", "Qupath_Nucleus_Hematoxylin_OD_mean",
                          "Qupath_Nucleus_Eosin_OD_mean", "Qupath_NucleustoCell_area_ratio", "Qupath_Cell_Area", "Qupath_Cell_Perimeter", "Qupath_Cell_Circularity",
                          "Qupath_Cell_Hematoxylin_OD_mean", "Qupath_Cell_Eosin_OD_mean", "Qupath_Cytoplasm_Hematoxylin_OD_mean", "Qupath_Cytoplasm_Eosin_OD_mean",
                          "gender", "age_sample",
                          "Splicing_mut", "DNA_chromatin_mut", "Cellcycle_mut", "Celldifferentiation_mut", "RAS_pathway_mut", mutations)
      
      survival_phenotypes <- c("AML", "sample_to_aml_or_last_date_date", "os_time", "os_event", "aza_duration_days_real_cutoff28days", "aza_stop_event")
      
      # Variable types
      ## Continuous variables
      continuous_phenotypes <- c("Stroma_proportion", "Lipid_droplet_proportion", "RBC_proportion", "WBC_proportion", "IPSS_score", "IPSSR_score",
                                 "cytogenetics_ipss_class_score", "cytogenetics_ipssr_class_score",
                                 "b-hb", "b-leuk", "b-neut", "b_trom", "bm_blast_p", "cytopenia_score_ipss", "cytopenia_score_ipssr",
                                 "Qupath_Nucleus_Area", "Qupath_Nucleus_Perimeter", "Qupath_Nucleus_Circularity", "Qupath_Nucleus_Eccentricity", "Qupath_Nucleus_Hematoxylin_OD_mean",
                                 "Qupath_Nucleus_Eosin_OD_mean", "Qupath_NucleustoCell_area_ratio", "Qupath_Cell_Area", "Qupath_Cell_Perimeter", "Qupath_Cell_Circularity",
                                 "Qupath_Cell_Hematoxylin_OD_mean", "Qupath_Cell_Eosin_OD_mean", "Qupath_Cytoplasm_Hematoxylin_OD_mean", "Qupath_Cytoplasm_Eosin_OD_mean",
                                 "age_sample")
      ## Categorical variables
      categorical_phenotypes <- c("Diagnosis", "Dysplasia", "IPSS_class", "IPSSR_class", "Complex_karyotype", "Del5q", "Karyotype", "Disease_status", "WHO",
                                  "secondary_MDS", "Ymonosomy", "Del11q", "Del12p", "Del20q", "Del7q", "Monosomy7", "Trisomy8", "Trisomy19", "i17q", "Ch3abn",
                                  "cytogenetics_ipss_class", "cytogenetics_ipssr_class",
                                  "gender",
                                  "Splicing_mut", "DNA_chromatin_mut", "Cellcycle_mut", "Celldifferentiation_mut", "RAS_pathway_mut", mutations)
      ### Remove variables with only one non-NA value
      categorical_phenotypes <- categorical_phenotypes[sapply(vgg16_feature_list_cluster[categorical_phenotypes], function(x) length(unique(x[!is.na(x)]))) > 1]
      
      
      ########### Intercluster analysis - All clusters in loop #######################################################################
      
      
      # Loop
      ## Continuous variables
      print("Loop over continuous variables")
      pvalue_df_cont <- data.frame()
      ignore_variables_df <- data.frame(matrix(ncol = 6, nrow = 0))
      ignore_variables_df1 <- data.frame(matrix(ncol = 6, nrow = 0))
      ignore_variables_df2 <- data.frame(matrix(ncol = 4, nrow = 0))
      ignore_variables_list <- list()
      for (clustering1 in clustering) {
        for (i in as.factor(unlist(unique(vgg16_feature_list_cluster[clustering1])))) {
          
          ## Make one cluster against the rest
          vgg16_feature_list_cluster$clusterX <- ifelse(is.na(vgg16_feature_list_cluster[clustering1]), NA,
                                                        ifelse(vgg16_feature_list_cluster[clustering1] == i, 1, 0))
          
          # Remove variables with all NA values
          ## Ignore variables with all NA in clusterX == 1
          ignore_variables <- lapply(vgg16_feature_list_cluster[continuous_phenotypes], function(x) max(x[vgg16_feature_list_cluster$clusterX==1], na.rm=TRUE))
          ignore_variables <- names(ignore_variables[ignore_variables == "-Inf"])
          ignore_variables_list <- c(ignore_variables_list, ignore_variables)
          ### Remove these variables
          continuous_phenotypes1 <- continuous_phenotypes[!continuous_phenotypes %in% ignore_variables]
          ## Ignore variables with all NA in clusterX == 0
          ignore_variables1 <- lapply(vgg16_feature_list_cluster[continuous_phenotypes1], function(x) max(x[vgg16_feature_list_cluster$clusterX==0], na.rm=TRUE))
          ignore_variables1 <- names(ignore_variables1[ignore_variables1 == "-Inf"])
          ignore_variables_list <- c(ignore_variables_list, ignore_variables1)
          ### Remove these variables
          continuous_phenotypes1 <- continuous_phenotypes1[!continuous_phenotypes1 %in% ignore_variables1]
          
          
          ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
          multiple_t_tests_p_value <- lapply(vgg16_feature_list_cluster[continuous_phenotypes1], function(x) wilcox.test(x ~ vgg16_feature_list_cluster$clusterX, na.rm=TRUE))
          # multiple_t_tests_p_value <- lapply(vgg16_feature_list_cluster["IPSSR_score"], function(x) wilcox.test(x ~ vgg16_feature_list_cluster$clusterX, na.rm=TRUE))
          ### P-values can be extracted from the result object
          pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
          ### Create a matrix and dataframe of the p-values
          pvalue_df1 <- pvalue %>% data.frame() %>% mutate(
            #### Add the p values to a new dataframe
            p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH"),
            #### Add also the t values, 95%CI to the same dataframe
            median1 = unlist(lapply(vgg16_feature_list_cluster[continuous_phenotypes1], function(x) median(x[vgg16_feature_list_cluster$clusterX==1], na.rm=TRUE))),
            median2 = unlist(lapply(vgg16_feature_list_cluster[continuous_phenotypes1], function(x) median(x[vgg16_feature_list_cluster$clusterX==0], na.rm=TRUE)))
          ) %>%
            ### Rownames to column
            rownames_to_column() %>%
            mutate(rowname = continuous_phenotypes1,
                   clusternumber = i) %>%
            rename(variable = rowname,
                   pvalue = p.value
            ) %>%
            arrange(pvalue)
          rownames(pvalue_df1) = NULL
          
          # Save data
          pvalue_df_cont <- rbind(pvalue_df_cont, pvalue_df1)
          
          # Add ignored variables as NA
          # ignore_variables2 <- c(ignore_variables, ignore_variables1)
          # ignore_variables_list
          for (h in ignore_variables_list) {
            ignore_variables_df <- data.frame("variable"  = h, "pvalue" = NA, "p_adjusted" = NA, "median1" = NA, "median2" = NA, "clusternumber" = i)
            ignore_variables_df1 <- rbind(ignore_variables_df1, ignore_variables_df)
          }
          
          colnames(ignore_variables_df1) <- colnames(pvalue_df_cont)
          pvalue_df_cont <- rbind(pvalue_df_cont, ignore_variables_df1)
          pvalue_df_cont <- pvalue_df_cont %>% filter(!is.na(variable)) %>% arrange(pvalue) %>% distinct(variable, clusternumber, .keep_all = TRUE)
        }
       
        ## Categorical variables
        print("Loop over categorical variables")
        pvalue_df_cat <- data.frame()
        ignore_variables_df2 <- data.frame(matrix(ncol = 4, nrow = 0))
        ignore_variables_df <- data.frame(matrix(ncol = 4, nrow = 0))
        ignore_variables_list <- list()
        pvalue_df_cat_observed_counts <- list()
        for (i in as.factor(unlist(unique(vgg16_feature_list_cluster[clustering1])))) {
          
          ## Make one cluster against the rest
          vgg16_feature_list_cluster$clusterX <- ifelse(is.na(vgg16_feature_list_cluster[clustering1]), NA,
                                                        ifelse(vgg16_feature_list_cluster[clustering1] == i, 1, 0))
          
          # Remove variables with all NA values
          ## Ignore variables with all NA in clusterX == 1
          ignore_variables <- lapply(vgg16_feature_list_cluster[categorical_phenotypes], function(x) max(as.numeric(as.character(x[vgg16_feature_list_cluster$clusterX==1])), na.rm=TRUE))
          ignore_variables <- names(ignore_variables[ignore_variables == "-Inf"])
          ignore_variables_list <- c(ignore_variables_list, ignore_variables)
          ### Remove these variables
          categorical_phenotypes1 <- categorical_phenotypes[!categorical_phenotypes %in% ignore_variables]
          ## Ignore variables with all NA in clusterX == 0
          ignore_variables1 <- lapply(vgg16_feature_list_cluster[categorical_phenotypes1], function(x) max(x[vgg16_feature_list_cluster$clusterX==0], na.rm=TRUE))
          ignore_variables1 <- names(ignore_variables1[ignore_variables1 == "-Inf"])
          ignore_variables_list <- c(ignore_variables_list, ignore_variables1)
          ### Remove these variables
          categorical_phenotypes1 <- categorical_phenotypes1[!categorical_phenotypes1 %in% ignore_variables1]
          
          
          ## Multiple comparison (two-sided, unpaired, Mann-Whitney U test)
          multiple_t_tests_p_value <- lapply(vgg16_feature_list_cluster[categorical_phenotypes1], function(x) chisq.test(vgg16_feature_list_cluster$clusterX, x))
          ### P-values can be extracted from the result object
          pvalue <- data.frame(p.value = sapply(multiple_t_tests_p_value, getElement, name = "p.value"))
          ### Observed counts
          observed_counts = lapply(multiple_t_tests_p_value, getElement, name = "observed")
          ### Create a matrix and dataframe of the p-values
          pvalue_df2 <- pvalue %>% data.frame() %>% mutate(
            #### Add the p values to a new dataframe
            p_adjusted = p.adjust(p = as.matrix(pvalue), method = "BH")
          ) %>%
            ### Rownames to column
            rownames_to_column() %>%
            mutate(rowname = categorical_phenotypes1,
                   clusternumber = i) %>%
            rename(variable = rowname,
                   pvalue = p.value
            ) %>%
            arrange(pvalue)
          rownames(pvalue_df2) = NULL
          
          
          # Add ignored variables as NA
          # ignore_variables2 <- c(ignore_variables, ignore_variables1)
          # ignore_variables_list
          for (h in ignore_variables_list) {
            ignore_variables_df <- data.frame("variable"  = h, "pvalue" = NA, "p_adjusted" = NA, "clusternumber" = i)
            ignore_variables_df2 <- rbind(ignore_variables_df2, ignore_variables_df)
          }
          
          pvalue_df_cat <- rbind(pvalue_df_cat, ignore_variables_df2)
          pvalue_df_cat <- pvalue_df_cat %>% filter(!is.na(variable)) %>% arrange(pvalue) %>% distinct(variable, clusternumber, .keep_all = TRUE)
          
          # Save data
          pvalue_df_cat <- rbind(pvalue_df_cat, pvalue_df2)
          pvalue_df_cat_observed_counts <- c(pvalue_df_cat_observed_counts, observed_counts)
        }
        
        # Export data
        print("Export data")
        ## Create directory
        if (file.exists(file.path("PCA_and_UMAP/Xception/Statistics"))){
          print("Export directory tree already exists")
        } else {
          print(paste0("Create ", "PCA_and_UMAP/Xception/Statistics"))
          dir.create(file.path("PCA_and_UMAP/Xception/Statistics"))
        }
        
        ## Export tables
        write_xlsx(pvalue_df_cont, paste0("PCA_and_UMAP/Xception/Statistics/Continuous_variable_", clustering1, "_statistics.xlsx"))
        write_xlsx(pvalue_df_cat, paste0("PCA_and_UMAP/Xception/Statistics/Categorical_variable_", clustering1, "_statistics_p.xlsx"))
        saveRDS(pvalue_df_cat_observed_counts, paste0("PCA_and_UMAP/Xception/Statistics/Categorical_variable_", clustering1, "_statistics_counts.rds"))
        
        
        ########### Balloon plot #######################################################################
        
        print("Generate balloon plot")
        df <- pvalue_df_cont[!is.na(pvalue_df_cont$pvalue),]
        # Modify data
        ## Rename variables
        unique(df$variable)
        df <- df %>% mutate(
          variable = gsub("Qupath_", "", variable),
          variable = gsub("b-leuk", "B-Leuk", variable),
          variable = gsub("b_trom", "B-PLT", variable),
          variable = gsub("b-neut", "B-Neut", variable),
          variable = gsub("b-leuk", "B-Leuk", variable),
          variable = gsub("b-hb", "B-Hb", variable),
          variable = gsub("cytogenetics_ipssr_class_score", "Cytogenetics score (IPSSR)", variable),
          variable = gsub("cytopenia_score_ipssr", "Cytopenia score (IPSSR)", variable),
          variable = gsub("cytogenetics_ipss_class_score", "Cytogenetics score (IPSS)", variable),
          variable = gsub("cytopenia_score_ipss", "Cytopenia score (IPSS)", variable),
          variable = gsub("age_sample", "Age", variable),
          variable = gsub("bm_blast_p", "BM Blasts %", variable),
          variable = gsub("_OD_", " ", variable),
          variable = gsub("_", " ", variable)
        )
        
        ## Calculate FC and -log10 P value
        df <- df %>%
          mutate(FC = median1 / median2) %>%
          mutate(
            neg_log10_p = -log10(pvalue),
            neg_log10_p_adj = -log10(p_adjusted),
            p_adjusted_cat = ifelse(p_adjusted<0.001, 0.001, 
                                    ifelse(p_adjusted<0.01, 0.01,
                                           ifelse(p_adjusted<0.05, 0.05,
                                                  ifelse(p_adjusted<0.1, 0.1, "ns")))),
            pvalue_cat = ifelse(pvalue<0.001, 0.001, 
                                ifelse(pvalue<0.01, 0.01,
                                       ifelse(pvalue<0.05, 0.05,
                                              ifelse(pvalue<0.1, 0.1, "ns")))),
            FC_log = log(FC) )
        
        
        # Balloonplot
        p <- ggballoonplot(df, x = "clusternumber", y = "variable",
                           fill = "FC_log",
                           size = "neg_log10_p_adj",
                           size.range = c(1, 10),
                           ggtheme = theme_bw()) +
          # scale_fill_viridis_c(option = "C") +
          gradient_fill(c("blue", "white", "red")) +
          guides(size = guide_legend(title="-log10(adjusted P value)"),
                 fill = guide_colorbar(title="Log10 fold change")) +
          font("xy.text", size = 10, color = "black", face="plain")
        # p
        
        # Export plot
        print("Export plot")
        png(paste0("PCA_and_UMAP/Xception/Statistics/Balloonplot_continuous_", clustering1, "_variables.png"), width = 6, height = 8, units = 'in', res = 300)
        print(p)
        dev.off()
        
      }
    }
  }
}
