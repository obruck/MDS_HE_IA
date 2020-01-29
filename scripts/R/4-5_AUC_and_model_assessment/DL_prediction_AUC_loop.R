# Created: OB, 01.12.2019

# This script assesses and explains elastic net models by visualizing prediction activations and other analyses
## 0. Defines parameters --> lines 0-74
## 1. Recalls models and the datasets (used in model training in a separate cluster) --> lines 75-539
## 2. Plots lambda, cve and beta coefficient values for the models --> lines 540-622
## 3. Generates ROCs and AUCs for the models at tile and TMA (WSI) level --> lines 623-1104
## 4. Generates representative TMA and tile images for each model --> lines 1105-1492
## 5. Collects AUC values for each model and plots the best values and parameteres into a heatmap --> lines 1488-1599
## 6. Generates a correlation matrix plot of the correlation between models and pixel- and cell-level image analysis data --> lines 1600-1698
## 7. Plots a dendrogram for the correlation matrix --> lines 1699-1710

# NB! Depending on your working directory paths, some paths and strings will need to be changed (especially in part 4, lines 1108-1492)


rm(list=ls())

# Load packages
library(tidyverse)
library(stringr)
library(readxl)
library(writexl)
library(dplyr)
library(glmnet)
library(biglasso)
library(pROC)
library(survAUC)
library(plotROC)
library(ROCR)
library(raster)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(ggdendro)
library(dendextend)


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


##################################### Feature matrix ##############################################################################################


for (w in w_list) {
  for (y in y_list) {
    # }}
    
    # Load Imagenet data
    if (y == "VGG16") {
      ## VGG16
      vgg16_feature_list1 <- readRDS("/data/processed_files/Feature_list/Feature_list_VGG16_MDS.rds")
      print(paste0("Load Imagenet data with parameters ", w, " and ", y))
    } else if (y == "Xception") {
      ## Xception
      vgg16_feature_list1 <- readRDS("/data/processed_files/Feature_list/Feature_list_Xception_MDS.rds")
      print(paste0("Load Imagenet data with parameters ", w, " and ", y))
    } else {
      print(paste0("Not enough parameters to load Imagenet data with parameters ", w, " and ", y))
    }
    
    
    
    ## Make a variable with the image spot name
    vgg16_feature_list1$name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles/", "", vgg16_feature_list1$vgg16_feature_list_name)
    vgg16_feature_list1$name <- gsub("\\/MDS_TMA", "MDS_TMA", vgg16_feature_list1$name)
    vgg16_feature_list1$name <- gsub("_modified[[:print:]]+", "", vgg16_feature_list1$name)
    vgg16_feature_list1 <- vgg16_feature_list1 %>% dplyr::select(name, everything())
    
    
    
    #################################### Clinical data ##############################################################################################
    
    
    # Load data
    print("Load clinical data")
    vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")
    
    
    # Keep only MDS, MDS/MPN and AA patients
    vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="MDS" |
                                                                                 vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="PreMDS",]
    
    
    # Merge clinical and feature data
    clinical <- vgg16_feature_list_cluster_clinical %>%
      data.frame() %>%
      dplyr::rename(name = Image_name)
    ## Define dataframe based on parameter w
    ## Define k = label to use in plots for the variable
    if (w == "Cmplx") {
      clinical <- clinical %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Complex karyotype")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "del5q") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = ch_five_q_del) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Ch5del")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "aml") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(event = os_event, time = sample_to_aml_time) %>%
        dplyr::select(name, event, time)
      assign("k", "Risk of AML")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "aml_log") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = AML_2y) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Risk of AML log")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "overall_survival") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(event = AML, time = sample_to_last_date_time) %>%
        dplyr::select(name, event, time)
      assign("k", "Risk of death")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "overall_survival_log") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = os_event_2y) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Risk of death log")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "aza") {
      vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="MDS" |
                                                                                   vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="PreMDS",]
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$Baseline_1ypreazastart_aza==1,]
      
      clinical <- vgg16_feature_list_cluster_clinical %>%
        data.frame() %>%
        dplyr::rename(name = Image_name) %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(event = aza_stop_event, time = aza_duration_days_real_cutoff28days) %>%
        dplyr::select(name, event, time)
      assign("k", "Aza treatment response")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "aza_log") {
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$Baseline_1ypreazastart_aza==1,]
      
      clinical <- vgg16_feature_list_cluster_clinical %>%
        data.frame() %>%
        dplyr::rename(name = Image_name) %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = aza_6mo_response) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Aza treatment response log")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "abnchr") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = diagnose_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Any chromosomal aberration")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "asxl1") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = ASXL1_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "ASXL1 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "Cellcycle_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Cellcycle_mut) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Mutation in Cell cycle pathway")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "Celldifferentiation_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Celldifferentiation_mut) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Mutation in cell differentiation pathway")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "ch3abn") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Ch3abn) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Chromosome 3 abnormality")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "del20q") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Del20q) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Del20q")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "del7q") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Del7q) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Del7q")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "dg") {
      vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="MDS" |
                                                                                   vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="PreMDS" |
                                                                                   vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample_HEproject=="MDS/MPN",]
      
      clinical <- vgg16_feature_list_cluster_clinical %>%
        data.frame() %>%
        dplyr::rename(name = Image_name) %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = diagnosis_abbr_sample_mds_mpdsmpn) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Diagnosis")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "DNA_chromatin_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = DNA_chromatin_mut) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Mutation in DNA chromatin pathway")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "gender") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = gender_num) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Gender")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "idh_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = IDH_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "IDH mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "monosomy7") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Monosomy7) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Monosomy 7")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "ras_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = RAS_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "NRAS/KRAS mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "RAS_pathway_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = RAS_pathway_mut) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Mutation in RAS pathway")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "runx1") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = RUNX1_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "RUNX1 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "secondary_MDS") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = secondary_MDS) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Secondary MDS")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "splicing_mut") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Splicing_mut) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Splicing mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "sf3b1") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = SF3B1_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "SF3B1 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "srsf2") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = SRSF2_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "SRSF2 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "stag2") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = STAG2_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "STAG2 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "tet2") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = TET2_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "TET2 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "tp53") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = TP53_kval) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "TP53 mutation")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "trisomy8") {
      clinical <- clinical %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = Trisomy8) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Trisomy 8")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "IPSSR") {
      clinical <- clinical %>%
        filter(Disease_status == "Dg") %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = ipssr_score) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "IPSSR score")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "age_dg") {
      clinical <- clinical %>%
        filter(Disease_status == "Dg") %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = age_dg) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "Age")
      print(paste0("Modify data using ", w, " data"))
    } else if (w == "cytopeniaipssr") {
      clinical <- clinical %>%
        filter(Disease_status == "Dg") %>%
        dplyr::select(-cmplx_ipss) %>%
        dplyr::rename(cmplx_ipss = cytopenia_score_ipssr) %>%
        dplyr::select(name, cmplx_ipss)
      assign("k", "IPSSR Cytopenia score")
      print(paste0("Modify data using ", w, " data"))
    } else {
      print("Error! Not enough parameters.")
    }
    
    
    # Left join clinical data to CNN feature matrix
    if (!w %in% w_list_cox & w %in% z_list) {
      hm1 <- clinical %>%
        left_join(vgg16_feature_list1, by = "name") %>%
        dplyr::select(name, vgg16_feature_list_name, cmplx_ipss, everything())
      print(paste0("Left join clinical data to CNN feature matrix using logistic or linear regression and ", w, " and z_list"))
    } else if (w %in% w_list_cox & w %in% z_list) {
      hm1 <- clinical %>%
        left_join(vgg16_feature_list1, by = "name") %>%
        dplyr::select(name, vgg16_feature_list_name, event, time, everything())
      print(paste0("Left join clinical data to CNN feature matrix using Cox regression and ", w, " and z_list"))  
    } else if (!w %in% w_list_cox & !w %in% z_list) {
      hm1 <- vgg16_feature_list1 %>%
        left_join(clinical, by = "name") %>%
        dplyr::select(name, vgg16_feature_list_name, cmplx_ipss, everything())
      print(paste0("Left join clinical data to CNN feature matrix using logistic or linear regression and ", w, " without z_list"))
    } else if (w %in% w_list_cox & !w %in% z_list) {
      hm1 <- vgg16_feature_list1 %>%
        left_join(clinical, by = "name") %>%
        dplyr::select(name, vgg16_feature_list_name, event, time, everything())
      print(paste0("Left join clinical data to CNN feature matrix using Cox regression and ", w, " without z_list"))
    } else {
      print("Error! No parameters defined for left joining clinical data to CNN feature matrix")
    }
    
    
    #################################### Split dataset ##############################################################################################
    
    
    # Glmnet
    print(paste0("Split dataset with data on ", w))
    
    if (!w %in% w_list_cox) {
      ## Dg classes
      hm1$cmplx_ipss <- as.numeric(hm1$cmplx_ipss)
      hm1 <- hm1[!is.na(hm1$cmplx_ipss),]
      hm1 <- hm1[!is.na(hm1$vgg16_feature_list_name),]
      
      
      ## Sample data to training and test set
      ### Sample by wsi number
      set.seed(101)
      n <- length(unique(hm1$name))
      if (w == "aza_log") {
        sample <- sample(unique(hm1$name), size = n * 1, replace = FALSE)
      } else {
        sample <- sample(unique(hm1$name), size = n * 0.67, replace = FALSE)
      }
      
      
      ## Train
      ### One sampling
      train_mds <- hm1 %>% filter(name %in% sample)
      
      
      ### X
      trainX_mds <- train_mds[,3:ncol(hm1)] %>% dplyr::select(-cmplx_ipss)
      trainX_mds <- sapply(trainX_mds, function(x) as.numeric(x))
      trainX_mds <- as.matrix(trainX_mds)
      ### Y
      if (w %in% w_list_log) {
        trainY_mds <- as.factor(as.matrix(train_mds["cmplx_ipss"]))   # For logistic regression
        print(paste0("Predict ", w, " with logistic regression and n=", length(unique(train_mds$name)), " unique patients in the training dataset"))
      } else if (w %in% w_list_lin) {
        trainY_mds <- as.matrix(train_mds["cmplx_ipss"])              # For linear regression
        print(paste0("Predict ", w, " with linear regression and n=", length(unique(train_mds$name)), " unique patients in the training dataset"))
      } else {
        print(paste0("ERROR? You are now predicting ", w, " with linear regression"))
      }
      
      
      
      
      ## Test
      ### One sampling
      test_mds <- hm1 %>% filter(!name %in% sample)
      
      
      ### X
      testX_mds <- test_mds[,3:ncol(hm1)] %>% dplyr::select(-cmplx_ipss)
      testX_mds <- sapply(testX_mds, function(x) as.numeric(x))
      testX_mds <- as.matrix(testX_mds)
      ### Y
      if (w %in% w_list_log) {
        testY_mds <- as.factor(as.matrix(test_mds["cmplx_ipss"]))   # For logistic regression
        print(paste0("Predict ", w, " with logistic regression and n=", length(unique(test_mds$name)), " unique patients in the test dataset"))
      } else if (w %in% w_list_lin) {
        testY_mds <- as.matrix(test_mds["cmplx_ipss"])              # For linear regression
        print(paste0("Predict ", w, " with linear regression and n=", length(unique(test_mds$name)), " unique patients in the test dataset"))
      } else {
        print(paste0("ERROR? You are now predicting ", w, " with linear regression"))
      }
      
      ## As big matrix
      trainX_mds.bm <- as.big.matrix(trainX_mds)
      testX_mds.bm <- as.big.matrix(testX_mds)
      
    } else if (w %in% w_list_cox) {
      
      # Remove NA
      hm1 <- hm1[!(is.na(hm1$vgg16_feature_list_name) | is.na(hm1$event) | is.na(hm1$time)),]
      
      # Recode patients alive/non-AML after 5 years as survivors/non-AML
      hm1 <- hm1 %>% mutate(
        event = ifelse(event == 1 & time > (5*365.25), 0, event)
      )
      
      ## Sample data to training and test set
      ### Sample by wsi number
      set.seed(101)
      n <- length(unique(hm1$name))
      sample <- sample(unique(hm1$name), size = n * 0.67, replace = FALSE)
      
      
      ## Train
      ### One sampling
      train_mds <- hm1 %>% filter(name %in% sample)
      
      ### X
      trainX_mds <- train_mds[,3:ncol(hm1)] %>% dplyr::select(-event, -time)
      trainX_mds <- sapply(trainX_mds, function(x) as.numeric(x))
      trainX_mds <- as.matrix(trainX_mds)
      ### Y
      if (w %in% w_list_cox) {
        trainY_mds <- as.matrix(train_mds %>% dplyr::select(event, time))     # For Cox regression
        print(paste0("Predict ", w, " with Cox regression and n=", length(unique(train_mds$name)), " unique patients in the training dataset"))
      } else {
        print(paste0("ERROR? You are now predicting ", w, " with ??? regression"))
      }
      
      
      ## Test
      ### One sampling
      test_mds <- hm1 %>% filter(!name %in% sample)
      
      
      ### X
      testX_mds <- test_mds[,3:ncol(hm1)] %>% dplyr::select(-event, -time)
      testX_mds <- sapply(testX_mds, function(x) as.numeric(x))
      testX_mds <- as.matrix(testX_mds)
      ### Y
      if (w %in% w_list_cox) {
        testY_mds <- as.matrix(test_mds %>% dplyr::select(event, time))     # For Cox regression
        print(paste0("Predict ", w, " with Cox regression and n=", length(unique(test_mds$name)), " unique patients in the test dataset"))
      } else {
        print(paste0("ERROR? You are now predicting ", w, " with ??? regression"))
      }
      
    } else {
      print("Error! No parameters defined for left joining clinical data to CNN feature matrix")
    }
    
    
    
    ########### Imagenet #########################################################################################
    
    # Load data
    test <- readRDS(paste0("/data/processed_files/Elastic_net_results/", w, "/", y, "/", snakecase::to_lower_camel_case(y), "_pred_", w, "_elastinet_biglasso_MDSonly.rds"))
    if (length(test) > 1) {
      print(paste0("Load elastic net models and define parameters using ", k, ", variable ", w, " and ", y))
      
      # Create export directory tree
      if (file.exists(file.path("/data/supervised_learning", w))){
        print("Export directory tree already exists")
      } else {
        print(paste0("Create ", "/data/supervised_learning/", w, " and copy directory tree"))
        dir.create(file.path("/data/supervised_learning", w))
        sapply(paste0("/data/supervised_learning/directory_tree/",
                      list.files("/data/supervised_learning/directory_tree/")),
               function(x) file.copy(from = x,
                                     to = paste0("/data/supervised_learning/", w, "/"),
                                     recursive=TRUE))
      }
      
    } else {
      print(paste0("Error! Elastic net models not loaded correctly with parameters ", k, ", variable ", w, " and ", y))
    }
    
    # Separate list into elements
    print("Separate list into elements")
    names(test) <- c("fit.bin.enet_min_cv", "fit.bin.enet_min", "fit.bin.enet_1se",
                     "fit.bin.lasso_min_cv", "fit.bin.lasso_min", "fit.bin.lasso_1se",
                     "fit.bin.ridge_min_cv", "fit.bin.ridge_min", "fit.bin.ridge_1se",
                     "pred_lasso_min", "pred_ridge_min", "pred_enet_min", "pred_lasso_1se",
                     "pred_ridge_1se", "pred_enet_1se")
    summary(test)
    enet_cv <- test[[1]]
    enet_min <- test[[2]]
    enet_1se <- test[[3]]
    lasso_cv <- test[[4]]
    lasso_min <- test[[5]]
    lasso_1se <- test[[6]]
    ridge_cv <- test[[7]]
    ridge_min <- test[[8]]
    ridge_1se <- test[[9]]
    pred_lasso_min <- test[[10]]
    pred_ridge_min <- test[[11]]
    pred_enet_min <- test[[12]]
    pred_lasso_1se <- test[[13]]
    pred_ridge_1se <- test[[14]]
    pred_enet_1se <- test[[15]]
    
    
    print("Generate elastic net lambda and cve plots")
    
    if (!w %in% w_list_cox) {
      ## Elastic net
      cvfit <- enet_cv
      ggsave(plot(cvfit), file=paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/lambda_enet.png", sep=""), width = 6, height = 5, units = "in", dpi = 300, limitsize = FALSE)
      print(paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/lambda_enet.png", sep=""))
      
      png(paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/beta_enet.png", sep=""), width = 6, height = 5, units = 'in', res = 300)
      plot(cvfit$fit) + abline(v = log(cvfit$lambda.min), col = 2, lty = 2)
      dev.off()
      print("Plot 1/3")
      
      ## Lasso
      cvfit <- lasso_cv
      ggsave(plot(cvfit), file=paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/lambda_lasso.png", sep=""), width = 6, height = 5, units = "in", dpi = 300, limitsize = FALSE)
      
      png(paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/beta_lasso.png", sep=""), width = 6, height = 5, units = 'in', res = 300)
      plot(cvfit$fit) + abline(v = log(cvfit$lambda.min), col = 2, lty = 2)
      dev.off()
      print("Plot 2/3")
      
      ## Ridge
      cvfit <- ridge_cv
      ggsave(plot(cvfit), file=paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/lambda_ridge.png", sep=""), width = 6, height = 5, units = "in", dpi = 300, limitsize = FALSE)
      
      png(paste0("/data/supervised_learning/", w, "/Glmnet/", y ,"/beta_ridge.png", sep=""), width = 6, height = 5, units = 'in', res = 300)
      plot(cvfit$fit) + abline(v = log(cvfit$lambda.min), col = 2, lty = 2)
      dev.off()
      print("Plot 3/3")
    } else {
      print(paste0("Produce these plots later for ", w_list_cox))
    }
    
    ########### Lasso vs ridge #########################################################################################
    
    
    # Run the analysis 6x for each parameters
    # Select 1. lasso, ridge or elasticnet and 2. lambda.min or lambda.1se
    
    fit_list <- list(lasso_min, lasso_1se, ridge_min, ridge_1se, enet_min, enet_1se)
    names(fit_list) <- c("lasso_min", "lasso_1se", "ridge_min", "ridge_1se", "enet_min", "enet_1se")
    penalty = NULL; alpha = NULL; lambda = NULL; j = NULL
    
    
    for (fit1 in names(fit_list)) {
      if (fit1 == "lasso_min") {
        if (!w %in% w_list_cox) {
          assign("penalty", "lasso"); assign("alpha", 1.0); assign("lambda", "lambda.min"); assign("j", "lasso_lambdamin"); fit = lasso_min
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "lasso"); assign("alpha", 1.0); assign("lambda", unique(lasso_cv$lambda.min)); assign("j", "lasso_lambdamin")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else if (fit1 == "lasso_1se") {
        if (!w %in% w_list_cox) {
          assign("penalty", "lasso"); assign("alpha", 1.0); assign("lambda", "lambda.1se"); assign("j", "lasso_lambda1se"); fit = lasso_1se
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "lasso"); assign("alpha", 1.0); assign("lambda", unique(lasso_cv$lambda.1se)); assign("j", "lasso_lambda1se")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else if (fit1 == "ridge_min") {
        if (!w %in% w_list_cox) {
          assign("penalty", "ridge"); assign("alpha", 0.0); assign("lambda", "lambda.min"); assign("j", "ridge_lambdamin"); fit = ridge_min
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "ridge"); assign("alpha", 0.0); assign("lambda", unique(ridge_cv$lambda.min)); assign("j", "ridge_lambdamin")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else if (fit1 == "ridge_1se") {
        if (!w %in% w_list_cox) {
          assign("penalty", "ridge"); assign("alpha", 0.0); assign("lambda", "lambda.1se"); assign("j", "ridge_lambda1se"); fit = ridge_1se
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "ridge"); assign("alpha", 0.0); assign("lambda", unique(ridge_cv$lambda.1se)); assign("j", "ridge_lambda1se")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else if (fit1 == "enet_1se") {
        if (!w %in% w_list_cox) {
          assign("penalty", "enet"); assign("alpha", 0.5); assign("lambda", "lambda.1se"); assign("j", "enet_lambda1se"); fit = enet_1se
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "enet"); assign("alpha", 0.5); assign("lambda", unique(enet_cv$lambda.1se)); assign("j", "enet_lambda1se")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else if (fit1 == "enet_min") {
        if (!w %in% w_list_cox) {
          assign("penalty", "enet"); assign("alpha", 0.5); assign("lambda", "lambda.min"); assign("j", "enet_lambdamin"); fit = enet_min
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        } else {
          assign("penalty", "enet"); assign("alpha", 0.5); assign("lambda", unique(enet_cv$lambda.min)); assign("j", "enet_lambdamin")
          print(paste0("Produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
        }
      } else {
        print(paste0("Error to produce data for each lambda and alpha parameters using ", "penalty = ", penalty, "; alpha = ", alpha, "; lambda = ", lambda, "; j = ", j))
      }
      
      
      # Perform AUC analyses
      ## Logistic and linear regression analyses
      if (!w %in% w_list_cox) {
        print("Perform AUC analyses for logistic and linear regression analyses")
        
        # Use significant coefficients to assign tile-wise predictions
        ## Train data
        pred <- predict(fit, trainX_mds.bm, type="response", penalty = penalty, alpha = alpha)
        ## Bind predictions to labels
        train_data <- cbind(as.matrix(pred), as.numeric(as.character(trainY_mds))) %>% as.data.frame()
        
        ## Tile-wise logistic regression model
        if (w %in% w_list_log) {
          model <- glm(V2 ~ V1, family=binomial(link="logit"), data=train_data)
        } else if (w %in% w_list_lin) {
          model <- glm(V2 ~ V1, family=gaussian(link = "identity"), data=train_data)
        } else {
          print("Error! with model family")
        }
        
        summary(model)
        roc_train_tile <- roc(train_data$V2, train_data$V1, plot = TRUE)
        roc_train_tile
        
        ### Combine predictions in the image-level
        train_data$imagename <- train_mds$name
        train_data$vgg16_feature_list_name <- train_mds$vgg16_feature_list_name
        train_data_summary <- train_data %>%
          ungroup() %>%
          mutate(mean_V1 = mean(V1, na.rm=TRUE),
                 quantile_V1 = quantile(V1, 0.90, na.rm=TRUE)) %>%
          group_by(mean_V1, imagename) %>%
          summarise(
            image_pred_mean = mean(V1),
            image_pred_prop = sum(V1 > mean_V1) / n()  # proportion over 75% quantile
          )
        train_data_summary <- train_data_summary %>% left_join(train_data, by = "imagename") %>% distinct(imagename, .keep_all = TRUE)
        
        ### Image-wise logistic regression
        if (w %in% w_list_log) {
          model_train <- glm(V2 ~ image_pred_mean, family=binomial(link="logit"), data=train_data_summary)
        } else if (w %in% w_list_lin) {
          model_train <- glm(V2 ~ image_pred_mean, family=gaussian(link = "identity"), data=train_data_summary)
        } else {
          print("Error! with model family")
        }
        
        ### ROC
        roc_train_wsi <- roc(train_data_summary$V2, train_data_summary$image_pred_mean, plot = TRUE)
        roc_train_wsi
        
        if (nrow(test_mds) > 0) {
          ## Test data
          pred <- predict(fit, testX_mds.bm, type="response", penalty=penalty, alpha = alpha)
          ## Bind predictions to labels
          test_data <- cbind(as.matrix(pred), as.numeric(as.character(testY_mds))) %>% as.data.frame()
          
          ## Tile-wise logistic regression model
          if (w %in% w_list_log) {
            model <- glm(V2 ~ V1, family=binomial(link="logit"), data=test_data)
          } else if (w %in% w_list_lin) {
            model <- glm(V2 ~ V1, family=gaussian(link = "identity"), data=test_data)
          } else {
            print("Error! with model family")
          }
          summary(model)
          roc_test_tile <- roc(test_data$V2, test_data$V1, plot = TRUE)
          roc_test_tile
          
          ### Combine predictions in the image-level
          test_data$imagename <- test_mds$name
          test_data$vgg16_feature_list_name <- test_mds$vgg16_feature_list_name
          test_data_summary <- test_data %>%
            ungroup() %>%
            mutate(mean_V1 = mean(V1, na.rm=TRUE),
                   # quantile_V1 = quantile(V1)[4]) %>%
                   quantile_V1 = quantile(V1, 0.90, na.rm=TRUE)) %>%
            group_by(mean_V1, imagename) %>%
            summarise(
              image_pred_mean = mean(V1, na.rm=TRUE),
              image_pred_prop = sum(V1 > mean_V1) / n()  # proportion over 75% quantile
            )
          test_data_summary <- test_data_summary %>% left_join(test_data, by = "imagename") %>% distinct(imagename, .keep_all = TRUE)
          ### Image-wise logistic regression
          if (w %in% w_list_log) {
            model_test <- glm(V2 ~ image_pred_mean, family=binomial(link="logit"), data=test_data_summary)
          } else if (w %in% w_list_lin) {
            model_test <- glm(V2 ~ image_pred_mean, family=gaussian(link = "identity"), data=test_data_summary)
          } else {
            print("Error! with model family")
          }
          ### ROC
          roc_test_wsi <- roc(test_data_summary$V2, test_data_summary$image_pred_mean, plot = TRUE)
          roc_test_wsi
          
          # Combine train and test data in tiles
          ## Tiles
          train_data$dataset <- "train"
          test_data$dataset <- "test"
          train_test_data <- rbind(train_data, test_data)
          ## WSI
          train_data_summary$dataset <- "train"
          test_data_summary$dataset <- "test"
          train_test_data_summary <- rbind(train_data_summary, test_data_summary)
          
          # Make a dataframe of the AUC values
          roc_df <- data.frame(as.numeric(roc_train_tile$auc), as.numeric(roc_test_tile$auc), as.numeric(roc_train_wsi$auc), as.numeric(roc_test_wsi$auc)) %>%
            rename(train_tile_auc = 1,
                   test_tile_auc = 2,
                   train_wsi_auc = 3,
                   test_wsi_auc = 4)
          rownames(roc_df) <- paste0(w, "_", y, "_", j)
          
        } else {
          print("Not enough data for a test dataset")
          
          # Combine train and test data in tiles
          ## Tiles
          train_data$dataset <- "train"
          train_test_data <- train_data
          ## WSI
          train_data_summary$dataset <- "train"
          train_test_data_summary <- train_data_summary
          
          # Make a dataframe of the AUC values
          roc_df <- data.frame(as.numeric(roc_train_tile$auc), as.numeric(roc_train_wsi$auc)) %>%
            rename(train_tile_auc = 1,
                   train_wsi_auc = 2)
          rownames(roc_df) <- paste0(w, "_", y, "_", j)
        }
        
        # Export
        write_xlsx(roc_df, paste0("/data/supervised_learning/", w, "/ROC/", y, "/ROC_", w, "_", y, "_", j, ".xlsx", sep=""))
        
      } else if (w %in% w_list_cox) {
        
        ## Cox regression analyses
        print("Perform AUC analyses for Cox regression analyses")
        
        # Use significant coefficients to assign tile-wise predictions
        ## Train data
        fit <- glmnet(trainX_mds, Surv(trainY_mds[,2], trainY_mds[,1]), alpha = alpha, family = "cox", lambda = lambda)
        pred <- predict(fit, trainX_mds, s = lambda, type = "response")
        
        ## Bind predictions to labels
        train_data <- cbind(as.matrix(pred), trainY_mds) %>% as.data.frame() %>% dplyr::rename(V1 = 1, V2 = event)
        
        ## Tile-wise logistic regression model
        model <- coxph(Surv(time, V2) ~ V1, data=train_data)
        summary(model)
        roc_train_tile <- roc(train_data$V2, train_data$V1, plot = TRUE)
        
        ### Combine predictions in the image-level
        train_data$imagename <- train_mds$name
        train_data$vgg16_feature_list_name <- train_mds$vgg16_feature_list_name
        train_data_summary <- train_data %>%
          ungroup() %>%
          mutate(mean_V1 = mean(V1, na.rm=TRUE),
                 quantile_V1 = quantile(V1, 0.90, na.rm=TRUE)) %>%
          group_by(mean_V1, imagename) %>%
          summarise(
            image_pred_mean = mean(V1),
            image_pred_prop = sum(V1 > mean_V1) / n()  # proportion over 75% quantile
          )
        train_data_summary <- train_data_summary %>% left_join(train_data, by = "imagename") %>% distinct(imagename, .keep_all = TRUE)
        ### Image-wise logistic regression
        model_train <- coxph(Surv(time, V2) ~ image_pred_mean, data=train_data_summary)
        ### ROC
        roc_train_wsi <- roc(train_data_summary$V2, train_data_summary$image_pred_mean, plot = TRUE)
        roc_train_wsi
        
        
        ## Test data
        if (nrow(test_mds) > 0) {
          pred <- predict(fit, testX_mds, s = lambda, type="response")
          ## Bind predictions to labels
          test_data <- cbind(as.matrix(pred), testY_mds) %>% as.data.frame() %>% dplyr::rename(V1 = 1, V2 = event)
          
          ## Tile-wise logistic regression model
          model <- coxph(Surv(time, V2) ~ V1, data=test_data)
          summary(model)
          roc_test_tile <- roc(test_data$V2, test_data$V1, plot = TRUE)
          ### Combine predictions in the image-level
          test_data$imagename <- test_mds$name
          test_data$vgg16_feature_list_name <- test_mds$vgg16_feature_list_name
          test_data_summary <- test_data %>%
            ungroup() %>%
            mutate(mean_V1 = mean(V1, na.rm=TRUE),
                   quantile_V1 = quantile(V1, 0.90, na.rm=TRUE)) %>%
            group_by(mean_V1, imagename) %>%
            summarise(
              image_pred_mean = mean(V1, na.rm=TRUE),
              image_pred_prop = sum(V1 > mean_V1) / n()  # proportion over 75% quantile
            )
          test_data_summary <- test_data_summary %>% left_join(test_data, by = "imagename") %>% distinct(imagename, .keep_all = TRUE)
          ### Image-wise logistic regression
          model_test <- coxph(Surv(time, V2) ~ image_pred_prop, data=test_data_summary)
          model_test
          ### ROC
          roc_test_wsi <- roc(test_data_summary$V2, test_data_summary$image_pred_mean, plot = TRUE)
          roc_test_wsi
          
          
          # Combine train and test data in tiles
          ## Tiles
          train_data$dataset <- "train"
          test_data$dataset <- "test"
          train_test_data <- rbind(train_data, test_data)
          ## WSI
          train_data_summary$dataset <- "train"
          test_data_summary$dataset <- "test"
          train_test_data_summary <- rbind(train_data_summary, test_data_summary)
          
          
          # Make a dataframe of the AUC values
          roc_df <- data.frame(as.numeric(roc_train_tile$auc), as.numeric(roc_test_tile$auc), as.numeric(roc_train_wsi$auc), as.numeric(roc_test_wsi$auc)) %>%
            rename(train_tile_auc = 1,
                   test_tile_auc = 2,
                   train_wsi_auc = 3,
                   test_wsi_auc = 4)
          rownames(roc_df) <- paste0(w, "_", y, "_", j)
          # Export
          write_xlsx(roc_df, paste0("/data/supervised_learning/", w, "/ROC/", y, "/ROC_", w, "_", y, "_", j, ".xlsx", sep=""))
          
        } else {
          print("No enough data for a test dataset")
          
          # Combine train and test data in tiles
          ## Tiles
          train_data$dataset <- "train"
          train_test_data <- train_data
          ## WSI
          train_data_summary$dataset <- "train"
          train_test_data_summary <- train_data_summary
          
          
          # Make a dataframe of the AUC values
          roc_df <- data.frame(as.numeric(roc_train_tile$auc), as.numeric(roc_train_wsi$auc)) %>%
            rename(train_tile_auc = 1,
                   train_wsi_auc = 2)
          rownames(roc_df) <- paste0(w, "_", y, "_", j)
          # Export
          write_xlsx(roc_df, paste0("/data/supervised_learning/", w, "/ROC/", y, "/ROC_", w, "_", y, "_", j, ".xlsx", sep=""))
          
        }
        
      } else {
        print("Error! Neither logistic, linear nor Cox regression analyses")
      }
      
      ########### ROC plots #########################################################################################
      
      
      theme_map <- function (base_size = 20, base_family = "") {
        theme_gray(base_size = base_size, base_family = base_family) %+replace%
          theme(legend.position=c(0.77, 0.155))
      }
      
      # Loop over both "tiles" and "wsi"
      images = c("tiles", "wsi")
      if (w != "aza_log") {
        for (i in images) {
          if (i == "tiles") {
            # Combine train and test data
            ## Tiles
            train_data$dataset <- "train"
            test_data$dataset <- "test"
            train_test_data <- rbind(train_data, test_data)
            
            ## ROC
            roc_train <- roc_train_tile
            roc_test <- roc_test_tile
            
            # Define predictions for Train and Test separately
            ## Tiles
            train_test_data$TrainProb <- ifelse(train_test_data$dataset=="train", train_test_data$V1, NA)
            train_test_data$TestProb <- ifelse(train_test_data$dataset=="test", train_test_data$V1, NA)
            
            # Yhdistä data
            ## Tile
            longtest <- melt_roc(train_test_data, "V2", c("TrainProb", "TestProb"))
            longtest <- longtest %>% filter(!is.na(D)) %>% rename(Dataset = name)
            
            # Export
            ## Logistic and Cox regression analyses only!
            if (!w %in% w_list_lin) {
              print("Generate AUC plots for logistic and Cox regression analyses")
              
              g1 <- ggplot(longtest, aes(d = D, m = M, color = Dataset)) +
                scale_color_brewer("Dataset", labels = c(paste("TestProb", round(roc_test$auc, 3), sep=" "), paste("TrainProb", round(roc_train$auc, 3), sep=" ")), palette = "Set1") +
                geom_roc(n.cuts = 0, show.legend = TRUE, size = 2) +
                ggtitle(k) +
                theme(plot.title = element_text(hjust = 0.5)) +
                style_roc(theme = theme_map, ylab = "Sensitivity", xlab = "1 - Specificity")  #theme = theme_map
              ggsave(g1, file=paste0("/data/supervised_learning/", w, "/ROC/", y, "/AUC_", i, "_", y, "_", j,".png"), width = 7, height = 5, limitsize = FALSE, units = 'in', dpi = 300)
            } else { print("No AUC plots for linear regression analyses") }
          } else {
            # Combine train and test data
            ## WSI
            train_data_summary$dataset <- "train"
            test_data_summary$dataset <- "test"
            train_test_data_summary <- rbind(train_data_summary, test_data_summary)
            
            ## ROC
            roc_train <- roc_train_wsi
            roc_test <- roc_test_wsi
            
            # Define predictions for Train and Test separately
            ## WSI
            train_test_data_summary$TrainProb <- ifelse(train_test_data_summary$dataset=="train", train_test_data_summary$image_pred_prop, NA)
            train_test_data_summary$TestProb <- ifelse(train_test_data_summary$dataset=="test", train_test_data_summary$image_pred_prop, NA)
            
            #yhdistä data
            ## WSI
            longtest <- melt_roc(train_test_data_summary, "V2", c("TrainProb", "TestProb"))
            longtest <- longtest %>% filter(!is.na(D.V2)) %>% rename(Dataset = name, D = D.V2) %>% dplyr::select(-D.V2.1)
            
            # Export
            ## Logistic and Cox regression analyses only!
            if (!w %in% w_list_lin) {
              print("Generate AUC plots for logistic and Cox regression analyses")
              
              g2 <- ggplot(longtest, aes(d = D, m = M, color = Dataset)) +
                scale_color_brewer("Dataset", labels = c(paste("TestProb", round(roc_test$auc, 3), sep=" "), paste("TrainProb", round(roc_train$auc, 3), sep=" ")), palette = "Set1") +
                geom_roc(n.cuts = 0, show.legend = TRUE, size = 2) +
                ggtitle(k) +
                theme(plot.title = element_text(hjust = 0.5)) +
                style_roc(theme = theme_map, ylab = "Sensitivity", xlab = "1 - Specificity")  #theme = theme_map
              ggsave(g2, file=paste0("/data/supervised_learning/", w, "/ROC/", y, "/AUC_", i, "_", y, "_", j,".png"), width = 7, height = 5, limitsize = FALSE, units = 'in', dpi = 300)
            } else { print("No AUC plots for linear regression analyses") }
          }
        }
      } else {
        for (i in images) {
          if (i == "tiles") {
            # Combine train and test data
            ## Tiles
            train_data$dataset <- "train"
            train_test_data <- train_data
            
            ## ROC
            roc_train <- roc_train_tile
            
            # Define predictions for Train and Test separately
            ## Tiles
            train_test_data$TrainProb <- ifelse(train_test_data$dataset=="train", train_test_data$V1, NA)
            
            # Yhdistä data
            ## Tile
            longtest <- melt_roc(train_test_data, "V2", c("TrainProb", "TrainProb"))
            longtest <- longtest %>% filter(!is.na(D)) %>% rename(Dataset = name) %>% filter(Dataset == "TrainProb")
            
            # Export
            ## Logistic and Cox regression analyses only!
            if (!w %in% w_list_lin) {
              print("Generate AUC plots for logistic and Cox regression analyses")
              
              g1 <- ggplot(longtest, aes(d = D, m = M, color = Dataset)) +
                scale_color_brewer("Dataset", labels = paste("TrainProb", round(roc_train$auc, 3), sep=" "), palette = "Set1") +
                geom_roc(n.cuts = 0, show.legend = TRUE, size = 2) +
                ggtitle(k) +
                theme(plot.title = element_text(hjust = 0.5)) +
                style_roc(theme = theme_map, ylab = "Sensitivity", xlab = "1 - Specificity")  #theme = theme_map
              ggsave(g1, file=paste0("/data/supervised_learning/", w, "/ROC/", y, "/AUC_", i, "_", y, "_", j,".png"), width = 7, height = 5, limitsize = FALSE, units = 'in', dpi = 300)
            } else { print("No AUC plots for linear regression analyses") }
          } else {
            # Combine train and test data
            ## WSI
            train_data_summary$dataset <- "train"
            train_test_data_summary <- train_data_summary
            
            ## ROC
            roc_train <- roc_train_wsi
            
            # Define predictions for Train and Test separately
            ## WSI
            train_test_data_summary$TrainProb <- ifelse(train_test_data_summary$dataset=="train", train_test_data_summary$image_pred_prop, NA)
            
            #yhdistä data
            ## WSI
            longtest <- melt_roc(train_test_data_summary, "V2", c("TrainProb", "TrainProb"))
            longtest <- longtest %>% filter(!is.na(D.V2)) %>% rename(Dataset = name, D = D.V2) %>% dplyr::select(-D.V2.1) %>% filter(Dataset == "TrainProb")
            
            # Export
            ## Logistic and Cox regression analyses only!
            if (!w %in% w_list_lin) {
              print("Generate AUC plots for logistic and Cox regression analyses")
              
              g2 <- ggplot(longtest, aes(d = D, m = M, color = Dataset)) +
                scale_color_brewer("Dataset", labels = paste("TrainProb", round(roc_train$auc, 3), sep=" "), palette = "Set1") +
                geom_roc(n.cuts = 0, show.legend = TRUE, size = 2) +
                ggtitle(k) +
                theme(plot.title = element_text(hjust = 0.5)) +
                style_roc(theme = theme_map, ylab = "Sensitivity", xlab = "1 - Specificity")  #theme = theme_map
              ggsave(g2, file=paste0("/data/supervised_learning/", w, "/ROC/", y, "/AUC_", i, "_", y, "_", j,".png"), width = 7, height = 5, limitsize = FALSE, units = 'in', dpi = 300)
            } else { print("No AUC plots for linear regression analyses") }
          }
        }
      }
      
      
      ########### Export raw data #########################################################################################
      
      
      ## Make list of results
      if (w != "aza_log") {
        imagenet_vgg16 <- list(pred, train_data, test_data, train_data_summary, test_data_summary, model, roc_train_tile, roc_test_tile, roc_train_wsi, roc_test_wsi)  # cut-off = 95% quantile of V1
      } else {
        imagenet_vgg16 <- list(pred, train_data, train_data_summary, model, roc_train_tile, roc_train_wsi)  # cut-off = 95% quantile of V1
      }
      
      # Save list
      saveRDS(imagenet_vgg16, paste0("/data/processed_files/Elastic_net_results/", w, "/models/imagenet_", y, "_", j, ".rds"))
      
    }
    
    
    ########### Make images from tiles #########################################################################################
    
    
    # Select best model and generate images with that
    
    # Define x and y coordinates of each tiles
    ## Combine train and test data
    if(nrow(test_mds) == 0) {
      print("No test dataset")
      test_mds <- train_mds[0,]
      test_data <- train_data[0,]
    } else {
      print("Test dataset ok")
    }
    mds <- rbind(train_mds, test_mds)
    ## X
    mds$x_coord <- gsub("[[:print:]]+jpg_", "", mds$vgg16_feature_list_name)
    mds$x_coord <- gsub("\\.jpg", "", mds$x_coord)
    mds$x_coord1 <- as.numeric(gsub("_[[:digit:]]*", "", mds$x_coord))
    ## Y
    mds$y_coord1 <- as.numeric(gsub("[[:digit:]]*_", "", mds$x_coord))
    ## Combination
    mds$x_coord <- paste(mds$x_coord1, mds$y_coord1, sep="_")
    
    # Combine predictions and coordinates
    train_test_data_pred <- mds %>% data.frame %>% dplyr::select(vgg16_feature_list_name, x_coord1, y_coord1, x_coord)
    ## Combine train and test data
    train_test_data <- rbind(train_data, test_data)
    train_test_data_pred <- train_test_data_pred %>% left_join(train_test_data, by = "vgg16_feature_list_name")
    
    ## Center and scale predictions to [0,1]
    train_test_data_pred$V1_scaled <- scale(train_test_data_pred$V1, center = TRUE, scale = TRUE)
    train_test_data_pred$V1_scaled <- as.numeric(train_test_data_pred$V1_scaled)
    
    ## Split by images into lists
    train_test_data_pred_by_image <- split(train_test_data_pred, f = train_test_data_pred$imagename)
    train_test_data_pred_by_image_1 <- train_test_data_pred_by_image[[2]]
    print(unique(train_test_data_pred_by_image_1[!is.na(train_test_data_pred_by_image_1$imagename),]$imagename))
    
    ## Create scaffold
    scaffold <- data.frame(matrix(nrow = 22*23, ncol = 2))
    scaffold$X1 <- rep(1:22, 23)
    scaffold$X2 <- rep(1:23, 22)
    scaffold <- scaffold %>% rename(x_coord1 = X1, y_coord1 = X2) %>% mutate(x_coord = paste(x_coord1, y_coord1, sep="_")) %>% arrange(x_coord1, y_coord1)
    
    
    # Merge predictions to scaffold
    train_test_data_pred_by_image_1 <- full_join(scaffold, train_test_data_pred_by_image_1) %>% arrange(x_coord1, y_coord1)
    
    
    # Plot
    ggplot(data=train_test_data_pred_by_image_1) + theme_bw() + coord_equal() +
      geom_raster(aes(x=y_coord1, y=-x_coord1, fill=V1_scaled)) +
      
      
      scale_fill_gradient2(low="#053061", mid="#e0e0e0", high = "#d6604d", na.value="white") +
      labs(fill = "Probability") +
      
      #General aesthetics
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            # legend.position="none",
            legend.key.size = unit(1, "cm"),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
      )
    
    
    ############ Generate plots for each image #########################################################################################
    
    
    for (i in 1:length(train_test_data_pred_by_image)) {
      train_test_data_pred_by_image_1 <- train_test_data_pred_by_image[[i]]
      print(unique(train_test_data_pred_by_image_1[!is.na(train_test_data_pred_by_image_1$imagename),]$imagename))
      
      ## Create scaffold
      scaffold <- data.frame(matrix(nrow = 22*23, ncol = 2))
      scaffold$X1 <- rep(1:22, 23)
      scaffold$X2 <- rep(1:23, 22)
      scaffold <- scaffold %>% rename(x_coord1 = X1, y_coord1 = X2) %>% mutate(x_coord = paste(x_coord1, y_coord1, sep="_")) %>% arrange(x_coord1, y_coord1)
      
      
      # Merge predictions to scaffold
      train_test_data_pred_by_image_1 <- full_join(scaffold, train_test_data_pred_by_image_1) %>% arrange(x_coord1, y_coord1)
      
      # Plot
      g <- ggplot(data=train_test_data_pred_by_image_1) + theme_bw() + coord_equal() +
        geom_raster(aes(x=y_coord1, y=-x_coord1, fill=V1_scaled)) +
        
        scale_fill_gradient2(low="#053061", mid="#e0e0e0", high = "#d6604d", na.value="white") +
        labs(fill = "Probability") +
        
        #General aesthetics
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.key.size = unit(1, "cm"),
              legend.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 12, face = "bold"),
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()
        )
      
      ggsave(g, file=paste0("/data/supervised_learning/", w, "/Prob/", y, "/", names(train_test_data_pred_by_image)[i],".png"), width = 10, height = 10, limitsize = FALSE, units = 'in', dpi = 300)
    }
    
    
    ############ Plots of top and bottom 100 tiles #########################################################################################
    
    
    if (!w %in% w_list_lin) {
      # Logistic and Cox regression
      
      ## Select top and bottom 100 tiles
      top_tiles <- train_test_data_pred %>% ungroup() %>% arrange(V1) %>% top_n(2000, V1) %>% sample_n(100)
      bottom_tiles <- train_test_data_pred %>% ungroup() %>% arrange(V1) %>% top_n(-2000, V1) %>% sample_n(100)
      
      ## Identify the folders
      to_folder_top <- paste0("/data/supervised_learning/", w, "/Tiles_examples/", y, "/Top/")
      to_folder_bottom <- paste0("/data/supervised_learning/", w, "/Tiles_examples/", y, "/Bottom/")
      ## Find the files that you want
      top_tiles$vgg16_feature_list_name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles//", "/data/images/3_preprocessed_tiles_grayscale/", top_tiles$vgg16_feature_list_name)
      file_list_top <- top_tiles$vgg16_feature_list_name
      bottom_tiles$vgg16_feature_list_name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles//", "/data/images/3_preprocessed_tiles_grayscale/", bottom_tiles$vgg16_feature_list_name)
      file_list_bottom <- bottom_tiles$vgg16_feature_list_name
      
      ## Copy the files to the new folder
      file.copy(from=file_list_top, to=to_folder_top,
                # overwrite = recursive, recursive = FALSE, 
                copy.mode = TRUE)
      file.copy(from=file_list_bottom, to=to_folder_bottom,
                copy.mode = TRUE)
      print(paste0("Plots of top and bottom 100 tiles for ", w))
    } else if (w %in% w_list_lin) {
      
      # Select top and bottom 100 tiles
      ipssr_tiles <- train_test_data_pred %>% arrange(V1) %>% mutate(V1 = as.numeric(V1)) %>%
        mutate(ipssr_group = Hmisc::cut2(x=V1, g=10)) %>% mutate(ipssr_group = cumsum(!duplicated(ipssr_group))) %>%
        group_by(ipssr_group) %>% sample_n(10)
      # Copy images
      to_folder_ipssr <- paste0("/data/supervised_learning/", w, "/Tiles_examples/", y, "/Top/")
      ## Find the files that you want
      ipssr_tiles$vgg16_feature_list_name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles//", "/data/images/3_preprocessed_tiles_grayscale/", ipssr_tiles$vgg16_feature_list_name)
      file_list_ipssr_tiles <- ipssr_tiles$vgg16_feature_list_name
      
      ## Copy the files to the new folder
      file.copy(from=file_list_ipssr_tiles, to=to_folder_ipssr,
                copy.mode = TRUE)
      
      
      ## Rename files in the new folder
      #### IPSSR
      ipssr_tiles$new_file_name <- gsub("\\/data\\/images\\/3_preprocessed_tiles_grayscale\\/",
                                        to_folder_ipssr, ipssr_tiles$vgg16_feature_list_name)
      ipssr_tiles$new_file_name1 <- gsub("\\.jpg", "", ipssr_tiles$new_file_name)
      files <- paste(ipssr_tiles$new_file_name1, "_Group", ipssr_tiles$ipssr_group, ".jpg", sep="")
      file.rename(from=ipssr_tiles$new_file_name,
                  to=files)
      
      
      # Plot linear regression
      for (dataseti in unique(train_test_data_summary$dataset)) {
        g <- train_test_data_summary %>% filter(dataset == dataseti) %>% ggplot(aes(x=image_pred_mean, y=V2)) +
          geom_point(col="steelblue", size=3) +
          geom_smooth(method = "lm", color="black") +
          labs(
            x=paste("Prediction probability of ", tolower(k), sep=""),
            y=gsub("i", "I", tolower(k))) +
          theme(axis.title.y = element_text(size=15,face="bold"),
                axis.title.x = element_text(size=15,face="bold")) +
          stat_cor(method = "pearson")
        ggsave(g, file=paste0("/data/supervised_learning/", w, "/Correlation/", y, "/Cor_", snakecase::to_snake_case(tolower(k)), "_pred_with_", dataseti, ".png", sep=""), width = 6, height = 6, limitsize = FALSE, units = 'in', dpi = 300)
      }
      print(paste0("Plots of top and bottom 100 tiles for ", w))
    } else {
      print(paste0("Error! Not enough data for plots of top and bottom 100 tiles for ", w))
    }
    
    
    ############ Plots of top and bottom 5 images #########################################################################################
    
    
    # Select top and bottom 5 wsis
    if (!w %in% w_list_lin) {
      top_wsi <- train_test_data_summary %>% ungroup() %>% arrange(-image_pred_mean) %>% filter(V2 == 1) %>% slice(1:5)
      bottom_wsi <- train_test_data_summary %>% ungroup() %>% arrange(image_pred_mean) %>% filter(V2 == 0) %>% slice(1:5)
      print(paste0("Select top and bottom 5 wsis with ", w))
    } else if (w %in% w_list_lin) {
      top_wsi <- train_test_data_summary %>% ungroup() %>% arrange(-image_pred_mean) %>% slice(1:5)
      bottom_wsi <- train_test_data_summary %>% ungroup() %>% arrange(image_pred_mean) %>% slice(1:5)
      print(paste0("Select top and bottom 5 wsis with ", w))
    } else {
      print(paste0("Error! Not enough data to select top and bottom 5 wsis with ", w))
    }
    
    # 1. Copy H&E images
    
    ## Identify the folders
    to_folder_top <- paste0("/data/supervised_learning/", w, "/Wsi_examples/", y,"/Top/")
    to_folder_bottom <- paste0("/data/supervised_learning/", w, "/Wsi_examples/", y,"/Bottom/")
    
    ## Find the files that you want
    top_wsi$vgg16_feature_list_name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles//", "/data/images/1_preprocessed_tma_color/", top_wsi$vgg16_feature_list_name)
    top_wsi$vgg16_feature_list_name <- gsub("\\.jpg[[:print:]]*", ".jpg", top_wsi$vgg16_feature_list_name)
    top_wsi$vgg16_feature_list_name <- gsub("compressed_size\\/", "compressed_size/Compressedsize_Compressed_", top_wsi$vgg16_feature_list_name)
    file_list_top <- top_wsi$vgg16_feature_list_name
    bottom_wsi$vgg16_feature_list_name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles//", "/data/images/1_preprocessed_tma_color/", bottom_wsi$vgg16_feature_list_name)
    bottom_wsi$vgg16_feature_list_name <- gsub("\\.jpg[[:print:]]*", ".jpg", bottom_wsi$vgg16_feature_list_name)
    bottom_wsi$vgg16_feature_list_name <- gsub("compressed_size\\/", "compressed_size/Compressedsize_Compressed_", bottom_wsi$vgg16_feature_list_name)
    file_list_bottom <- bottom_wsi$vgg16_feature_list_name
    
    ## Copy the files to the new folder
    file.copy(from=file_list_top, to=to_folder_top,
              copy.mode = TRUE)
    file.copy(from=file_list_bottom, to=to_folder_bottom,
              copy.mode = TRUE)
    
    
    ## Rename files in the new folder
    ### Remove "Compressedsize_Compressed_"
    #### TOP
    files <- list.files(paste0("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/"))
    file.rename(from=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/", files, sep=""),
                to=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/", sub(pattern="Compressedsize_Compressed_", replacement="", files), sep=""))
    #### BOTTOM
    files <- list.files(paste0("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/"))
    file.rename(from=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/", files, sep=""),
                to=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/", sub(pattern="Compressedsize_Compressed_", replacement="", files), sep=""))
    ### Remove "_modified"
    #### TOP
    files <- list.files(paste0("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/"))
    file.rename(from=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/", files, sep=""),
                to=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Top/", sub(pattern="_modified", replacement="", files), sep=""))
    #### BOTTOM
    files <- list.files(paste0("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/"))
    file.rename(from=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/", files, sep=""),
                to=paste("/data/supervised_learning/", w, "/Wsi_examples/", y, "/Bottom/", sub(pattern="_modified", replacement="", files), sep=""))
    print(paste0("Plots of top and bottom 5 images for ", w, ". H&E images copied."))
    
    
    # 2.Copy predictions
    
    ## Find the files that you want
    top_wsi$prediction_name <- gsub("2_MDS_HE[[:print:]]*", paste0("supervised_learning/", w, "/Prob/", y, "/"), top_wsi$vgg16_feature_list_name)
    top_wsi$prediction_name <- paste(top_wsi$prediction_name, top_wsi$imagename, ".png", sep="")
    file_list_top <- top_wsi$prediction_name
    bottom_wsi$prediction_name <- gsub("2_MDS_HE[[:print:]]*", paste0("supervised_learning/", w, "/Prob/", y, "/"), bottom_wsi$vgg16_feature_list_name)
    bottom_wsi$prediction_name <- paste(bottom_wsi$prediction_name, bottom_wsi$imagename, ".png", sep="")
    file_list_bottom <- bottom_wsi$prediction_name
    
    
    ## Copy the files to the new folder
    file.copy(from=file_list_top, to=to_folder_top,
              copy.mode = TRUE)
    file.copy(from=file_list_bottom, to=to_folder_bottom,
              copy.mode = TRUE)
    print(paste0("Plots of top and bottom 5 images for ", w, ". Prediction images copied."))
    
    
    # 3. Copy pixel classification images
    ## Find the files that you want
    file_list_bottom_px <- gsub("supervised_learning/[[:print:]]*/Prob/VGG16/", "4_weka_images/Compressedsize_Compressed_", file_list_bottom)
    file_list_bottom_px <- gsub("\\.png", "_modified.jpg", file_list_bottom_px)
    file_list_top_px <- gsub("supervised_learning/[[:print:]]*/Prob/VGG16/", "4_weka_images/Compressedsize_Compressed_", file_list_top)
    file_list_top_px <- gsub("\\.png", "_modified.jpg", file_list_top_px)
    
    
    ## Copy the files to the new folder
    file.copy(from=file_list_top_px, to=to_folder_top,
              copy.mode = TRUE)
    file.copy(from=file_list_bottom_px, to=to_folder_bottom,
              copy.mode = TRUE)
    
    ## Rename files
    ### TOP
    files <- list.files(to_folder_top)[str_detect(list.files(to_folder_top), "Compressed")]
    file.rename(from=paste(to_folder_top, files, sep=""),
                to=paste(to_folder_top, sub(pattern="\\.jpg", replacement="_px.jpg", files), sep=""))
    files <- list.files(to_folder_top)[str_detect(list.files(to_folder_top), "Compressed")]
    file.rename(from=paste(to_folder_top, files, sep=""),
                to=paste(to_folder_top, sub(pattern="Compressedsize_Compressed_", replacement="", files), sep=""))
    ### BOTTOM
    files <- list.files(to_folder_bottom)[str_detect(list.files(to_folder_bottom), "Compressed")]
    file.rename(from=paste(to_folder_bottom, files, sep=""),
                to=paste(to_folder_bottom, sub(pattern="\\.jpg", replacement="_px.jpg", files), sep=""))
    files <- list.files(to_folder_bottom)[str_detect(list.files(to_folder_bottom), "Compressed")]
    file.rename(from=paste(to_folder_bottom, files, sep=""),
                to=paste(to_folder_bottom, sub(pattern="Compressedsize_Compressed_", replacement="", files), sep=""))
    print(paste0("Plots of top and bottom 5 images for ", w, ". Pixel classification images copied."))
    
    
    ############ Correlate with pixel-classification data #######################################################################
    
    
    # Load data
    vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")
    vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>% dplyr::select(
      patient_id, Image_name, "stroma_proportion", "tumor_proportion", "rbc_proportion", "lipid_droplet_proportion_of_mask",
      ch_five_q_del, bm_blast_p, ipssr_score, cytopenia_score_ipss, cytopenia_score_ipssr, 
      contains("anemia_score_ipss"), contains("neutropenia_score_ipss"), contains("trombosytopenia_score_ipss"),
      contains("Qupath")) %>%
      rename("imagename" = "Image_name",
             "anemia_score_ipss" = "anemia_score_ipss =IF($M60<80; 1.5; IF($M60<100; 1; IF($M60>99.9; 0; \"\")))",
             "neutropenia_score_ipss" = "neutropenia_score_ipss =IF($O60<0.8; 0.5; IF($O60>0.799; 0; \"\"))",
             "trombosytopenia_score_ipss" = "trombosytopenia_score_ipss =IF($P60<50; 1; IF($P60<100; 0.5; IF($P60>99.9; 0; \"\")))",
             lipid_droplet = lipid_droplet_proportion_of_mask,
             red_blood_cells = rbc_proportion,
             white_blood_cells = tumor_proportion,
             stroma = stroma_proportion
      )
    
    # Merge with prediction data
    train_test_data_summary_clinical <- train_test_data_summary %>% left_join(vgg16_feature_list_cluster_clinical, by = "imagename")
    
    
    # Make separate
    cor_matrix <- train_test_data_summary_clinical %>%
      ungroup() %>%
      dplyr::select(image_pred_mean, white_blood_cells, red_blood_cells, lipid_droplet, stroma, contains("Qupath"))
    colnames(cor_matrix) <- sapply(colnames(cor_matrix), function(x) gsub("Qupath_", "", x))
    ## Export
    write_xlsx(cor_matrix, paste0("/data/supervised_learning/", w, "/Correlation/", y, "/Cor_", w, ".xlsx", sep=""))
    ## As.matrix
    cor_matrix <- cor_matrix %>% as.matrix()
    
    
    # Remove Qupath from colnames
    colnames(train_test_data_summary_clinical) <- sapply(colnames(train_test_data_summary_clinical), function(x) gsub("Qupath_", "", x))
    
    
    
    
    # Plot
    print(paste0("Generate correlation plots for ", w, " using pixel classification data."))
    pixel_classes1 = c("stroma", "red_blood_cells", "white_blood_cells", "lipid_droplet")
    pixel_classes2 = gsub("Qupath_", "", colnames(vgg16_feature_list_cluster_clinical)[str_detect(colnames(vgg16_feature_list_cluster_clinical), "Qupath")])
    for (pixel_class in pixel_classes1) {
      g <- ggplot(train_test_data_summary_clinical, aes_string(x="image_pred_mean", y=pixel_class)) +
        geom_point(col="steelblue", size=3) +
        geom_smooth(method = "lm", color="black") +
        labs(
          x=paste("Prediction probability of ", tolower(k), sep=""),
          y=gsub("_", " ", paste("Proportion of ", pixel_class, sep=""))) +
        theme(axis.title.y = element_text(size=15,face="bold"),
              axis.title.x = element_text(size=15,face="bold")) +
        stat_cor(method = "spearman")
      ggsave(g, file=paste0("/data/supervised_learning/", w, "/Correlation/", y, "/Cor_", snakecase::to_snake_case(tolower(k)), "_pred_with_", pixel_class, ".png", sep=""), width = 6, height = 6, limitsize = FALSE, units = 'in', dpi = 300)
    }
    
    print(paste0("Generate correlation plots for ", w, " using QuPath analytical data."))
    for (pixel_class in pixel_classes2) {
      g <- ggplot(train_test_data_summary_clinical, aes_string(x="image_pred_mean", y=pixel_class)) +
        geom_point(col="steelblue", size=3) +
        geom_smooth(method = "lm", color="black") +
        labs(
          x=paste("Prediction probability of ", tolower(k), sep=""),
          y=gsub("_", " ", pixel_class)) +
        theme(axis.title.y = element_text(size=15,face="bold"),
              axis.title.x = element_text(size=15,face="bold")) +
        stat_cor(method = "spearman")
      ggsave(g, file=paste0("/data/supervised_learning/", w, "/Correlation/", y, "/Cor_", snakecase::to_snake_case(tolower(k)), "_pred_with_", pixel_class, ".png", sep=""), width = 6, height = 6, limitsize = FALSE, units = 'in', dpi = 300)
    }
    
    rm(trainX_mds.bm); rm(testX_mds.bm); rm(trainX_mds); rm(testX_mds); rm(train_data); rm(test_data); rm(train_test_data); rm(pred); rm(fit)
  }
}


######################################## COMBINE AUCS #########################################################################################################


# Define alpha and lambda parameters
j_list = c("lasso_lambdamin", "lasso_lambda1se", "ridge_lambdamin", "ridge_lambda1se", "enet_lambdamin", "enet_lambda1se")
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
  data.frame()

# Make matrix for corrplot and add rownames
## With parameter information
colnames(auc_files_best)[1:4] <- c("Train AUC (Tiles)", "Test AUC (Tiles)", "Train AUC (WSI)", "Test AUC (WSI)")
rownames(auc_files_best) <- paste("AUC", auc_files_best$auc_files_param_y, auc_files_best$auc_files_param_alpha, auc_files_best$auc_files_param_lambda, sep="_")
auc_files_best_mat_w <- as.matrix(auc_files_best)[,1:4]
auc_files_best_mat_w <- apply(auc_files_best_mat_w, 2, as.numeric)
rownames(auc_files_best_mat_w) <- rownames(auc_files_best)
## Without parameter information
rownames(auc_files_best) <- snakecase::to_upper_lower_case(auc_files_best$auc_files_param_y)
auc_files_best_mat_wo <- as.matrix(auc_files_best)[,1:4]
auc_files_best_mat_wo <- apply(auc_files_best_mat_wo, 2, as.numeric)
rownames(auc_files_best_mat_wo) <- rownames(auc_files_best)
rownames(auc_files_best_mat_wo) <- gsub("CMPLX", "ComplexKaryotype", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("DG", "MDS_vs_MDS/MPN", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("TP53", "TP53mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("STAG2", "STAG2mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("ASXL1", "ASXL1mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("RUNX1", "RUNX1mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("TET2", "TET2mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("SF3b1", "SF3B1mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("SRSF2", "SRSF2mut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("RASpathwayMUT", "RASPATHWAYmut", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("SECONDARYmds", "SECONDARY_MDS", rownames(auc_files_best_mat_wo))
rownames(auc_files_best_mat_wo) <- gsub("DNAchromatinMUT", "DNACHROMATINmut", rownames(auc_files_best_mat_wo))


# AUC matrix plot
## With parameter information
png("/data/supervised_learning/AUC_matrix_with_parameters.png", width = 15, height = 15, units = 'in', res = 300)
corrplot(auc_files_best_mat_w,
         order = "original",
         cl.lim=c(0.5,1), # label bar range
         col = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),  #(200)
         addCoef.col = "black", # Cor values
         tl.cex = 1.5, # label font size
         number.cex = 1, # cor font size
         is.corr = FALSE, # Non-correlation matrix
         cl.ratio = 0.5, # width of colorlabel
         method="color",  # squares
         cl.align.text = "c",
         outline = TRUE,  # outline around squares
         tl.col = "black") # color of outline around squares
dev.off()
## Without parameter information
png("/data/supervised_learning/AUC_matrix_without_parameters.png", width = 15, height = 15, units = 'in', res = 300)
corrplot(auc_files_best_mat_wo,
         order = "original",
         cl.lim=c(0.5,1), # label bar range
         col = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),  #(200)
         addCoef.col = "black", # Cor values
         tl.cex = 1.5, # label font size
         number.cex = 1, # cor font size
         is.corr = FALSE, # Non-correlation matrix
         cl.ratio = 0.5, # width of colorlabel
         method="color",  # squares
         cl.align.text = "c",
         outline = TRUE,  # outline around squares
         tl.col = "black") # color of outline around squares
dev.off()



############ Correlate with pixel-classification data #######################################################################


# Correlation matrix for loop
## Load data
w1 = list(); w2 = list()
w_list <- c(w_list, "aza_log", w_list_lin)
for (w in w_list) {
  for (y in y_list) {
  i <- paste0("/data/supervised_learning/", w, "/Correlation/", y, "/Cor_", w, ".xlsx", sep="")
  w1 = c(w1, i)
  }
  w2 = c(w2, w1)
  w2 = unique(w2)
}
w3 = lapply(w2, read_xlsx)


## Make correlation matrix
## Correlation matrix cor values
cor_matrix <- lapply(w3, function(x) rcorr(as.matrix(x), type="spearman"))
## Correlation matrix cor values
cor_matrix_r <- lapply(cor_matrix, function(x) rbind(x$r[1,2:ncol(x$r)]))
cor_matrix_r <- do.call(rbind.data.frame, cor_matrix_r)
## Define rownames and convert to matrix
w2 <- gsub("[[:print:]]*Correlation\\/", "", w2); w2 <- gsub("\\.xlsx", "", w2); w2 <- gsub("\\/", "_", w2)
### Select only values in AUC matrix + for linear regression Xception values as Xception performed better than VGG16 for these
### Identify first whether to use Xception or VGG16
w2_tmp <- c(rownames(auc_files_best_mat_w), "age_dg", "cytopeniaipssr", "IPSSR")
w2_tmp <- w2_tmp[order(w2_tmp)]
w2_tmp <- ifelse(str_detect(w2_tmp, "VGG16"), paste0("VGG16_Cor_", w2_tmp), paste0("Xception_Cor_", w2_tmp))
w2_tmp <- gsub("_AUC", "", w2_tmp); w2_tmp <- gsub("_Xception_[[:print:]]*", "", w2_tmp); w2_tmp <- gsub("_VGG16_[[:print:]]*", "", w2_tmp)
w2_tmp <- c(w2_tmp, "Xception_Cor_aza_log")
### Cor matrix
cor_matrix_r <- cor_matrix_r %>%
  mutate(rowname = w2) %>%
  filter(rowname %in% w2_tmp) %>%
  column_to_rownames(var = "rowname") %>%
  as.matrix()
rownames(cor_matrix_r) <- gsub("[[:print:]]*Cor_", "", rownames(cor_matrix_r)) 
rownames(cor_matrix_r) <- gsub("_log", "", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- snakecase::to_upper_lower_case(rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("CMPLX", "ComplexKaryotype", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("DG", "MDS_vs_MDS/MPN", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("TP53", "TP53mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("STAG2", "STAG2mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("ASXL1", "ASXL1mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("RUNX1", "RUNX1mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("TET2", "TET2mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("SF3b1", "SF3B1mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("SRSF2", "SRSF2mut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("RASpathwayMUT", "RASPATHWAYmut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("DNAchromatinMUT", "DNACHROMATINmut", rownames(cor_matrix_r))
rownames(cor_matrix_r) <- gsub("SECONDARYmds", "SECONDARY_MDS", rownames(cor_matrix_r))
colnames(cor_matrix_r) <- gsub("white_blood_cells", "WBC", colnames(cor_matrix_r))
colnames(cor_matrix_r) <- gsub("red_blood_cells", "RBC", colnames(cor_matrix_r))
colnames(cor_matrix_r) <- gsub("lipid_droplet", "Lipid_droplet", colnames(cor_matrix_r))
colnames(cor_matrix_r) <- gsub("stroma", "Stroma", colnames(cor_matrix_r))
## Correlation matrix p values
cor_matrix_p <- lapply(cor_matrix, function(x) rbind(x$P[1,2:ncol(x$P)]))
cor_matrix_p <- do.call(rbind.data.frame, cor_matrix_p)
## Define rownames and convert to matrix
cor_matrix_p1 <- cor_matrix_p %>%
  dplyr::select(-c(contains("caliper"), contains("Perimeter"), contains("Circularity"))) %>%
  mutate(rowname = w2) %>%
  filter(rowname %in% w2_tmp) %>%
  column_to_rownames(var = "rowname") %>%
  sapply(function(x) p.adjust(x, method = "BH")) %>%
  as.matrix()


# Cluster rows by distance
cor_matrix_r1 <- cor_matrix_r %>% as.data.frame() %>% dplyr::select(-c(contains("caliper"), contains("Perimeter"), contains("Circularity"))) %>% as.matrix()
cor_matrix_r_distance <- dist(cor_matrix_r1, method = "euclidean")   # distance matrix
cor_matrix_r_hclust <- hclust(cor_matrix_r_distance, method="ward.D2") 
cor_matrix_r1 <- cor_matrix_r1[cor_matrix_r_hclust$order,]
cor_matrix_p1 <- cor_matrix_p1[cor_matrix_r_hclust$order,]

# Correlation matrix plot
png("/data/supervised_learning/Correlation_matrix2.png", width = 11, height = 9, units = 'in', res = 300)
corrplot(cor_matrix_r1, p.mat=cor_matrix_p1,
         is.corr = FALSE,
         sig.level = c(0.001, 0.01, 0.05),    # p-value sign limits
         insig = "label_sig",    # p-value sign
         col = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),  #(200)
         hclust.method = "ward.D2",  # clustering method
         pch=8,
         cl.length = 5,    # number of ticks in color key bar
         cl.lim = c(-1, 1),  # The limits(x1, x2)in the colorlabel
         pch.cex = 1.1,   # p-value sign size
         tl.srt = 90,  # label angle
         method="color",  # squares
         cl.align.text = "c",
         outline = TRUE,  # outline around squares
         tl.col = "black") # color of outline around squares
dev.off()


############ Dendrogram for corplot #######################################################################


png("/data/supervised_learning/Correlation_matrix1_dendrogram_plot2.png", width = 2, height = 4, units = 'in', res = 300)
dend <- as.dendrogram(cor_matrix_r_hclust) %>%
  set("branches_lwd", 1) %>% 
  set("labels_cex", 0.3) %>%   # plot2
  rotate((nrow(cor_matrix_r_hclust$merge)+1):1) %>%
  # set("labels_cex", 0.01) %>%    # plot
  plot(horiz=TRUE, axes=FALSE)
dev.off()

