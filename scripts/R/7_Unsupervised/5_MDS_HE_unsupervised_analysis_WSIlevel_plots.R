# Created: 31.03.2019, OB

# Plot UMAPs and PCAs

rm(list=ls())

# Load libraries
library(tidyverse)
library(Rphenograph)
library(uwot)
library(readxl)
library(RColorBrewer)
library(ggimage)
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
image_size1 <- c("wsi")


############### Load feature matrix data (pre-extracted from keras model) ###############################################################################


# Load feature matrix where the tile-level features have been averaged at the WSI-level (or average them for Xception)
for (dg_list in dg_list1) {
  for (image_size in image_size1) {
    for (dg_phase in dg_phase1) {
      # }}}
      set.seed(123)
      
      
      # Change working directory
      if (image_size == "wsi") {
        setwd(paste0("/data/results/unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"))
      } else {
        setwd(paste0("/data/results/unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"))
      }
      
      getwd()
      
      # Import data from first run
      vgg16_feature_list_cluster1 <- readRDS(paste0("Unsupervised_learning_", dg_list, ".rds"))
      
      
      ########### Plot UMAPs in batch mode - Set parameters ##############################################################################
      
      
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
      
      
      
      ########### Plot UMAPs in batch mode - Export plots #######################################################################
      
      
      # Plot UMAPs and PCAs in batch mode
      
      # Plotting parameters
      ## Directory
      i = paste(dg_list, dg_phase, sep="_")
      directory = dg_list
      variable_names <- c("Diagnosis", "Dysplasia", "IPSS_class", "IPSSR_class", "Complex_karyotype", "Del5q", "Karyotype", "Disease_status", "WHO",
                          "Stroma_proportion", "Lipid_droplet_proportion", "RBC_proportion", "WBC_proportion", "IPSS_score", "IPSSR_score")
      print(directory)
      
      
      ## X
      x <- c("PC1", "X1")
      ## Clustering
      clustering <- c("kcluster_pca", "phenograph_cluster_pca", "kcluster_umap", "phenograph_cluster_umap")
      clustering_kmeans <- c("kcluster_pca", "kcluster_umap")
      clustering_pheno <- c("phenograph_cluster_pca", "phenograph_cluster_umap")
      ## Variables
      continuous_phenotypes <- c("Stroma_proportion", "Lipid_droplet_proportion", "RBC_proportion", "WBC_proportion", "IPSS_score", "IPSSR_score")
      categorical_phenotypes <- c("Diagnosis", "Dysplasia", "IPSS_class", "IPSSR_class", "Complex_karyotype", "Del5q", "Karyotype", "Disease_status", "WHO")
      
      
      ## PCA + Kmeans on PCs
      for (i in variable_names) {
        if (i %in% continuous_phenotypes) {
          for (k in x) {
            if (k == "PC1") {
              for (m in clustering) {
                if (m %in% clustering_kmeans) {
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="PC2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_distiller(palette = "RdBu") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "K-means cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/", gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("PC1", "_onPCAcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                } else {  # m in clustering_pheno
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="PC2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_distiller(palette = "RdBu") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "Phenograph cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/", gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("PC1", "_onPCAcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                }
              }
            } else {   # k in UMAP
              for (m in clustering) {
                if (m %in% clustering_kmeans) {
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="X2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_distiller(palette = "RdBu") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "K-means cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("X1", "_onUMAPcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                } else {  # m in clustering_pheno
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="X2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_distiller(palette = "RdBu") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "Phenograph cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("X1", "_onUMAPcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                }
              }
            }
          }
        } else if (i %in% categorical_phenotypes) {
          for (k in x) {
            if (k == "PC1") {
              for (m in clustering) {
                if (m %in% clustering_kmeans) {
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="PC2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_brewer(palette = "Set1") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "K-means cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("PC1", "_onPCAcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                } else {  # m in clustering_pheno
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="PC2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_brewer(palette = "Set1") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "Phenograph cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("PC1", "_onPCAcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                }
              }
            } else {  # k in UMAP
              for (m in clustering) {
                if (m %in% clustering_kmeans) {
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="X2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_brewer(palette = "Set1") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "K-means cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("X1", "_onUMAPcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                } else {  # m in clustering_pheno
                  print(paste("Plot for ", i, " ", m, sep=""))
                  vgg16_feature_list_cluster$state <- vgg16_feature_list_cluster[,i]
                  g <- ggplot(vgg16_feature_list_cluster, aes_string(x=k, y="X2", shape = m)) +
                    geom_point(aes(color = state), size = 3) +
                    scale_shape_manual(values=c(8, 13, 3, 15, 16, 18, 22, 21, 24)) +
                    scale_color_brewer(palette = "Set1") +
                    stat_ellipse(aes_string(group = m)) +
                    labs(shape = "Phenograph cluster", color = gsub("_", " ", i))
                  ggsave(g, file=paste0("PCA_and_UMAP/", y, "/", i, "/",gsub("_", "", snakecase::to_snake_case(tolower(i))), "_Clustering_based_on_", snakecase::to_snake_case(tolower(m)), gsub("X1", "_onUMAPcoordinates", k), ".png", sep=""), width = 10, height = 7, limitsize = FALSE, units = 'in', dpi = 300)
                }
              }
            }
          }
        } else {
          print("Nothing to plot")
        }
      }
      
      
      ############# Find most representable WSI per cluster ##########################################################################################################
      
      
      # Visualize any point of the images
      vgg16_feature_list_cluster_identify <- vgg16_feature_list_cluster %>%
        dplyr::select(Image_name, everything())
      
      
      # Select point closest to cluster 1 median X1 and median X2
      ## Define median X1 and median X2
      example1 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==1) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example1$X1median, example1$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example1 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      
      
      # Select point closest to cluster 2 median X1 and median X2
      ## Define median X1 and median X2
      example2 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==2) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example2$X1median, example2$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example2 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      
      # Select point closest to cluster 3 median X1 and median X2
      ## Define median X1 and median X2
      example3 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==3) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example3$X1median, example3$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example3 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      
      # Select point closest to cluster 4 median X1 and median X2
      ## Define median X1 and median X2
      example4 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==4) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example4$X1median, example4$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example4 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      
      # Select point closest to cluster 5 median X1 and median X2
      ## Define median X1 and median X2
      example5 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==5) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example5$X1median, example5$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example5 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      
      # Select point closest to cluster 6 median X1 and median X2
      ## Define median X1 and median X2
      example6 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==6) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example6$X1median, example6$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example6 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      # Select point closest to cluster 7 median X1 and median X2
      ## Define median X1 and median X2
      example7 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==7) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example7$X1median, example7$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example7 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      # Select point closest to cluster 8 median X1 and median X2
      ## Define median X1 and median X2
      example8 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==8) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example8$X1median, example8$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example8 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      # Select point closest to cluster 9 median X1 and median X2
      ## Define median X1 and median X2
      example9 <- vgg16_feature_list_cluster_identify %>%
        dplyr::filter(phenograph_cluster_umap==9) %>%
        summarise(X1median = round(median(X1, na.rm=TRUE), 3),
                  X2median = round(median(X2, na.rm=TRUE), 3))
      ## Select point closest to (X1median, X2median)
      ## X1 and X2 coordinates
      p1 <- cbind(vgg16_feature_list_cluster_identify$X1, vgg16_feature_list_cluster_identify$X2)
      ## Coordinate to compare with
      p2 <- cbind(example9$X1median, example9$X2median)
      p2 <- p2[rep(seq_len(nrow(p2)), nrow(vgg16_feature_list_cluster)), ]
      ## Distance
      dm <- pointDistance(p1, p2, lonlat=FALSE, allpairs=FALSE)
      dm <- dm %>% as.data.frame() %>%
        rename(distance = ".") %>%
        mutate(name = vgg16_feature_list_cluster_identify$Image_name)
      ## Identify image
      example9 <- dm[dm$distance == min(dm$distance),] %>%
        dplyr::select(name)
      
      # print(c(example1, example2, example3, example4, example5, example6, example7, example8, example9))
      print(c(example1, example2, example3, example4, example5, example6, example7, example8))
      wsi_example <- list(example1, example2, example3, example4, example5, example6, example7)
      write.csv(wsi_example, "wsi_example.csv")
      
    }
  }
}