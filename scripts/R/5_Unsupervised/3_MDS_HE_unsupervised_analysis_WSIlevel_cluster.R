# Created: 31.03.2019, OB

# MDS HE clustering data in unsupervised way

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
image_size1 <- c("wsi", "tile")


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
      
      
      # Sample only 10 tiles per WSI
      if (image_size == "wsi") {
        print("wsi")
      } else {
        print("tile")
        vgg16_feature_list_cluster1 <- vgg16_feature_list_cluster1 %>%
          ungroup() %>%
          group_by(Image_name2) %>%
          sample_n(10) %>%
          ungroup()
      }
      
      
      ################## PCA ########################################################################
      
      
      # PCA
      print("PCA")
      vgg16_feature_list_cluster <- vgg16_feature_list_cluster1[,2:2049]
      pca <- prcomp(vgg16_feature_list_cluster,
                    center = FALSE,
                    scale = FALSE)
      
      
      ## Cumulative variance
      eigs <- pca$sdev^2
      test <- as.data.frame(t(rbind(
        SD = sqrt(eigs),
        Proportion = eigs/sum(eigs),
        Cumulative = cumsum(eigs)/sum(eigs))))
      
      ### Plot
      png("PCA_eigenvalues/Xception/PCA_eigenvalues_proportion.png", width = 6, height = 5, units = 'in', res = 300)
      g <- ggplot(test, aes(y=Cumulative, x=Proportion)) +
        geom_point(size = 3) +
        labs(y = "Cumulative Eigenvalues", x = "Eigenvalue of individual PC") +
        ylim(0,1)
      print(g)
      dev.off()
      
      test$rownumber <- as.data.frame(row(test)) %>% dplyr::select(V1) %>% unlist() %>% as.numeric()
      png("PCA_eigenvalues/Xception/PCA_eigenvalues_PC.png", width = 6, height = 5, units = 'in', res = 300)
      g <- ggplot(test, aes(y=Cumulative, x=rownumber)) +
        geom_point(size = 3) +
        labs(y = "Cumulative Eigenvalues", x = "Principal component") +
        ylim(0,1)
      print(g)
      dev.off()
      
      
      # Feature matrix from PCA
      vgg16_feature_list_cluster_pca <- as.data.frame(pca$x)
      
      # Filter PCs explaining <1% of the feature matrix variance
      ## Define Eigenvalues
      pca$eigs <- pca$sdev^2
      ## PC variation explanation = Eigenvalue/Total Eigenvalue
      pca$prop <- pca$eigs/sum(pca$eigs)
      ## Filter
      vgg16_feature_list_cluster_pca <- as.data.frame(rbind(vgg16_feature_list_cluster_pca, pca$prop))
      vgg16_feature_list_cluster_pca <- vgg16_feature_list_cluster_pca[, vgg16_feature_list_cluster_pca[nrow(vgg16_feature_list_cluster_pca),]>0.01]
      ## Remove Eigenvalues
      vgg16_feature_list_cluster_pca <- vgg16_feature_list_cluster_pca[1:(nrow(vgg16_feature_list_cluster_pca)-1),]
      
      
      ################### UMAP ########################################################################
      
      
      # Umap on PCs
      print("UMAP")
      umap_panel1 <- vgg16_feature_list_cluster1[,2:2049]
      
      # Compute umap
      umap_coord <- umap_panel1 %>%
        as.matrix() %>%
        umap()
      
      
      ################## KMeans ######################################################################################################
      
      
      # Optimizing k value for K-means
      
      print("Optimizing k value for K-means")
      
      # For PCA data
      
      # Elbow method
      png("K-means_kvalue_optimization/Xception/Screeplot_Kmeans_on_PCs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(vgg16_feature_list_cluster_pca, kmeans, method = "wss") +
        #geom_vline(linetype = 2, xintercept = 4) +
        labs(subtitle = "Elbow method")
      print(g)
      dev.off()
      
      # Silhouette method
      png("K-means_kvalue_optimization/Xception/Silhouette_Kmeans_on_PCs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(vgg16_feature_list_cluster_pca, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")
      print(g)
      dev.off()
      
      # Gap statistic
      # nboot = 50 to keep the function speedy.
      # recommended value: nboot= 500 for your analysis.
      # Use verbose = FALSE to hide computing progression.
      set.seed(123)
      if (image_size == "wsi") { nboots = 200 } else {nboots = 200}
      png("K-means_kvalue_optimization/Xception/Gapstat_Kmeans_on_PCs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(vgg16_feature_list_cluster_pca, kmeans, nstart = 10, k.max = 8, method = "gap_stat", nboot = nboots) +
        labs(subtitle = "Gap statistic method")
      print(g)
      dev.off()
      
      # 30 indices for choosing the best number of clusters
      png("K-means_kvalue_optimization/Xception/NbClust_Kmeans_on_PCs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- NbClust(data = vgg16_feature_list_cluster_pca, diss = NULL, distance = "euclidean",
                   min.nc = 3, max.nc = 8, method = "kmeans")  # method = “ward.D”, “ward.D2”, “single”, “complete”, “average”, “kmeans” and more.
      print(g)
      dev.off()
      center_pca = tail(names(sort(table(as.numeric(g$Best.nc)))), 1)
      
      
      
      ## For UMAP data
      
      # Elbow method
      png("K-means_kvalue_optimization/Xception/Screeplot_Kmeans_on_UMIs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(umap_coord, kmeans, method = "wss") +
        #geom_vline(linetype = 2, xintercept = 4) +
        labs(subtitle = "Elbow method")
      print(g)
      dev.off()
      
      # Silhouette method
      png("K-means_kvalue_optimization/Xception/Silhouette_Kmeans_on_UMIs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(umap_coord, kmeans, method = "silhouette")+
        labs(subtitle = "Silhouette method")
      print(g)
      dev.off()
      
      
      # Gap statistic
      # nboot = 50 to keep the function speedy.
      # recommended value: nboot= 500 for your analysis.
      # Use verbose = FALSE to hide computing progression.
      png("K-means_kvalue_optimization/Xception/Gapstat_Kmeans_on_UMIs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- fviz_nbclust(umap_coord, kmeans, nstart = 10, k.max = 8, method = "gap_stat", nboot = nboots) +
        labs(subtitle = "Gap statistic method")
      print(g)
      dev.off()
      
      
      # 30 indices for choosing the best number of clusters
      png("K-means_kvalue_optimization/Xception/NbClust_Kmeans_on_UMIs.png", width = 6, height = 5, units = 'in', res = 300)
      g <- NbClust(data = umap_coord, diss = NULL, distance = "euclidean",
                   min.nc = 3, max.nc = 8, method = "kmeans")  # method = “ward.D”, “ward.D2”, “single”, “complete”, “average”, “kmeans” and more.
      print(g)
      dev.off()
      
      center_umap = tail(names(sort(table(as.numeric(g$Best.nc)))), 1)
      
      
      ###################### Cluster ####################################################################
      
      
      print("Cluster data")
      
      # 1. Cluster PCs
      
      # A. K-means
      set.seed(123)
      kclust_pca <- kmeans(vgg16_feature_list_cluster_pca, centers=as.numeric(center_pca), nstart = 10)
      # B. Phenograph
      set.seed(123)
      phenoclust_pca <- Rphenograph(vgg16_feature_list_cluster_pca)
      
      
      # 2. Cluster UMIs
      
      # A. K-means
      ## K-means clustering solution
      set.seed(123)
      kclust_umap <- kmeans(umap_coord, centers=as.numeric(center_umap), nstart = 10)
      # B. Phenograph
      set.seed(123)
      phenoclust_umap <- Rphenograph(umap_coord)
      
      
      # 2. Cluster UMIs
      
      # A. K-means
      ## K-means clustering solution
      set.seed(123)
      kclust_fm <- kmeans(vgg16_feature_list_cluster1[,2:2049], centers=as.numeric(center_pca), nstart = 10)
      # B. Phenograph
      set.seed(123)
      phenoclust_fm <- Rphenograph(vgg16_feature_list_cluster1[,2:2049])
      
      
      ################### Export data for plots ###########################################################################
      
      
      print("Export data")
      
      # Add new clusters
      vgg16_feature_list_cluster1$kcluster_pca <- as.factor(kclust_pca$cluster)
      vgg16_feature_list_cluster1$kcluster_umap <- as.factor(kclust_umap$cluster)
      vgg16_feature_list_cluster1$kclust_fm <- as.factor(kclust_fm$cluster)
      vgg16_feature_list_cluster1$phenoclust_cluster_fm <- factor(membership(phenoclust_fm[[2]]))
      vgg16_feature_list_cluster1$phenoclust_cluster_pca <- factor(membership(phenoclust_pca[[2]]))
      vgg16_feature_list_cluster1$phenoclust_cluster_umap <- factor(membership(phenoclust_umap[[2]]))
      
      # Export data
      if (image_size == "wsi") {  
        saveRDS(vgg16_feature_list_cluster1, paste0("Unsupervised_learning_", dg_list, ".rds"))
        print("wsi")
      } else { 
        saveRDS(vgg16_feature_list_cluster1, paste0("Unsupervised_learning_", dg_list, "sample10perwsi.rds"))
        print("tile")
      }
    }
  }
}