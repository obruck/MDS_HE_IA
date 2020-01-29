# Created: 31.03.2019, OB

# Unsupervised learning of MDS HE TMA images

rm(list=ls())
library(keras)
library(tidyverse)
library(dplyr)
library(Rphenograph)
library(uwot)
library(readxl)
library(readr)
library(stringr)
library(RColorBrewer)
library(ggimage)
library(raster)
library(factoextra)
library(NbClust)


############### Parameters ###############################################################################


# Define which CNN network dataset to use
y_list <- "Xception"

# Define which diseases to include
dg_list1 <- c("MDS", "MDS_MDSMPN", "MDS_MDSMPN_CO", "MDS_MDSMPN_AA_CO", "MDS_CO")

# Define which disease stage to analyze
dg_phase1 <- c("dg", "dg_fu")

# Define WSI or tile
image_size <- c("wsi")

# Define z
z = "First round"


############### Load feature matrix data (pre-extracted from keras model) ###############################################################################


# Load feature matrix where the tile-level features have been averaged at the WSI-level (or average them for Xception)
for (dg_list in dg_list1) {
  for (dg_phase in dg_phase1) {
    if (!z == "First round over!" & y == "Xception") {
      
      
      ################# Load data ##########################################################################
      
      
      print("Load Xception data")
      vgg16_feature_list <- readRDS("/data/processed_files/Feature_list/Feature_list_Xception_MDS.rds")
      vgg16_feature_list$name <- gsub("\\/csc\\/mustjoki\\/imageanalysis\\/analysisR\\/input\\/mds_he_tiles\\/\\/", "", vgg16_feature_list$vgg16_feature_list_name)
      vgg16_feature_list$name <- gsub("_modified[[:print:]]*", "", vgg16_feature_list$name)
      vgg16_feature_list <- vgg16_feature_list %>%
        dplyr::select(-vgg16_feature_list_name) %>%
        dplyr::select(name, everything())
      
      ## Combine duplicate spots from the same control patient
      vgg16_feature_list_cluster_clinical <- read_xlsx("/Users/oscarbruck/OneDrive - University of Helsinki/Tutkimus/Projekteja/MDS/Analysis/HE/Data/Clinicaldata_forimages.xlsx") %>%
        dplyr::select(patient_id, Image_name, Image_name2) %>%
        mutate(name = Image_name)
      vgg16_feature_list <- vgg16_feature_list %>%
        left_join(vgg16_feature_list_cluster_clinical) %>%
        mutate(name = ifelse(str_detect(name, "CO_TMA"), Image_name2, name)) %>%
        dplyr::select(-patient_id, -Image_name2) %>%
        dplyr::select(name, everything())
      vgg16_feature_list <- vgg16_feature_list %>% dplyr::select(-Image_name)
      
      ## Aggregate by the image spot name
      vgg16_feature_list <- aggregate(.~name, vgg16_feature_list, mean)
      
      # Select numeric data
      vgg16_feature_list_cluster <- vgg16_feature_list %>% dplyr::select(-name)
      
      # Center and scale
      vgg16_feature_list_cluster_scaled <- scale(vgg16_feature_list_cluster, center = TRUE, scale = TRUE) 
    } else if (z == "First round over!") {
      print("Skip reloading data")
    } else {
      print("Error somewhere with z and y parameters")
    }
    
    
    ################# Feature engineering ##########################################################################
    
    
    set.seed(123)
    
    # Reload data faster
    vgg16_feature_list_cluster <- as.data.frame(vgg16_feature_list_cluster_scaled)
    
    
    # Save vgg16 feature list
    vgg16_feature_list_cluster1 <- vgg16_feature_list_cluster
    vgg16_feature_list_cluster <- vgg16_feature_list_cluster1
    
    
    ################## Load clinical data #########################################################################
    
    
    # Load clinical data
    vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")
    
    
    # Filter patients
    ## By diagnosis
    if (dg_list == "MDS") {
      print(paste0("Filter ", dg_list, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(diagnosis_abbr_sample_HEproject == "MDS" |
                 diagnosis_abbr_sample_HEproject=="PreMDS") %>%
        filter(!is.na(patient_id))
    } else if (dg_list == "MDS_MDSMPN") {
      print(paste0("Filter ", dg_list, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(diagnosis_abbr_sample_HEproject == "MDS" |
                 diagnosis_abbr_sample_HEproject=="PreMDS" |
                 diagnosis_abbr_sample_HEproject=="MDS/MPN") %>%
        filter(!is.na(patient_id))
    } else if (dg_list == "MDS_MDSMPN_CO") {
      print(paste0("Filter ", dg_list, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(diagnosis_abbr_sample_HEproject == "MDS" |
                 diagnosis_abbr_sample_HEproject=="PreMDS" |
                 diagnosis_abbr_sample_HEproject=="MDS/MPN" |
                 diagnosis_abbr_sample_HEproject=="Healthy") %>%
        filter(!is.na(patient_id))
    } else if (dg_list == "MDS_MDSMPN_AA_CO") {
      print(paste0("Filter ", dg_list, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(diagnosis_abbr_sample_HEproject == "MDS" |
                 diagnosis_abbr_sample_HEproject=="PreMDS" |
                 diagnosis_abbr_sample_HEproject=="MDS/MPN" |
                 diagnosis_abbr_sample_HEproject=="AA" |
                 diagnosis_abbr_sample_HEproject=="Healthy") %>%
        filter(!is.na(patient_id))
    } else if (dg_list == "MDS_CO") {
      print(paste0("Filter ", dg_list, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(diagnosis_abbr_sample_HEproject == "MDS" |
                 diagnosis_abbr_sample_HEproject=="PreMDS" |
                 diagnosis_abbr_sample_HEproject=="Healthy") %>%
        filter(!is.na(patient_id))
    } else {
      print("Error! Define which patients to analyze")
    }
    
    
    ## By disease stage
    vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>% mutate(Disease_status = factor(Disease_status, levels = c("Dg", "Pre", "FU")))
    if (dg_phase == "dg_fu") {
      print(paste0("Filter ", dg_phase, " phase samples"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>% filter(Disease_status == "Dg" | Disease_status == "Pre" | Disease_status == "FU")
    } else if (dg_phase == "dg") {
      print(paste0("Filter ", dg_phase, " patients"))
      vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical %>%
        filter(Disease_status == "Dg" | Disease_status == "Pre") %>% 
        arrange(patient_id, Disease_status) %>%
        group_by(patient_id) %>% 
        slice(1)
    } else {
      print("Error! Define which disease stage data to analyze")
    }
    
    
    # Modify data
    vgg16_feature_list_cluster <- vgg16_feature_list_cluster %>% mutate(name = vgg16_feature_list$name)
    
    # Merge clinical data
    vgg16_feature_list_cluster <- merge(vgg16_feature_list_cluster, vgg16_feature_list_cluster_clinical, all.x=TRUE, by="Image_name")
    
    
    ################## Create directory tree #########################################################################
    
    
    # Load data
    # Create export directory tree
    
    if (file.exists(file.path("/data/results/unsupervised_learning", paste(image_size, y, sep="_"), dg_phase, dg_list))){
      print("Export directory tree already exists")
    } else {
      print(paste0("Create ", "/data/results/unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, " and copy directory tree"))
      dir.create(file.path("/data/results/unsupervised_learning", paste(image_size, y, sep="_")))
      dir.create(file.path("/data/results/unsupervised_learning", paste(image_size, y, sep="_"), dg_phase))
      dir.create(file.path("/data/results/unsupervised_learning", paste(image_size, y, sep="_"), dg_phase, dg_list))
      sapply(paste0("/data/results/unsupervised_learning/directory_tree/",
                    list.files("/data/results/unsupervised_learning/directory_tree/")),
             function(x) file.copy(from = x,
                                   to = paste0("/data/results/unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"),
                                   recursive=TRUE
             ))
    }
    
    
    # Change working directory
    setwd(paste0("/data/results/unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"))
    getwd()
    
    
    ########################## Export data ##################################################################
    
    
    # Export and import data for plots
    saveRDS(vgg16_feature_list_cluster, paste0("Unsupervised_learning_", dg_list, ".rds"))
    
    z = "First round over!"
  }
}