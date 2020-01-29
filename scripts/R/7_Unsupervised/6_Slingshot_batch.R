# Created Nov 2019, OB
# Script for plotting slingshots ontop of umaps

rm(list=ls())

# Load packages
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(Rphenograph)
library(uwot)
library(SingleCellExperiment)


################ Parameters ##################################################################################################


# Define which CNN network dataset to use
y_list <- "Xception"

# Define which diseases to include
dg_list1 <- c("MDS", "MDS_MDSMPN", "MDS_MDSMPN_CO", "MDS_MDSMPN_AA_CO", "MDS_CO")

# Define which disease stage to analyze
dg_phase1 <- c("dg", "dg_fu")

# Define WSI or tile
# image_size <- c("wsi")
image_size <- c("tile")

# Define clustering methods to use
clustering <- c("kcluster_umap", "phenograph_cluster_umap")


################ Load data ##################################################################################################


# Load data
## WSI
for (dg_list in dg_list1) {
  for (y in y_list) {
    for (dg_phase in dg_phase1) {
      for (m in clustering) {
        set.seed(123)
        # }}}}

        # Change working directory
        if (image_size == "wsi") {
          setwd(paste0("/data/results/Unsupervised_learning/", paste(image_size, y, sep="_"), "/", y, "/", dg_phase, "/", dg_list, "/"))
        } else {
          setwd(paste0("/data/results/Unsupervised_learning/", paste(image_size, y, sep="_"), "/", dg_phase, "/", dg_list, "/"))
        }
        
        print(getwd())
        
        
        # Import data
        if (image_size == "wsi") {
          vgg16_feature_list_cluster <- readRDS(paste0("Unsupervised_learning_", dg_list, ".rds"))
          vgg16_feature_list_cluster <- vgg16_feature_list_cluster[order(vgg16_feature_list_cluster[,m]),]
          print("image_size == wsi")
        } else if (image_size == "tile") {
          vgg16_feature_list_cluster <- readRDS(paste0("Unsupervised_learning_", dg_list, "sample10perwsi.rds")) %>%
            as.data.frame() %>%
            dplyr::rename(phenograph_cluster_umap = phenoclust_cluster_umap)
          vgg16_feature_list_cluster <- vgg16_feature_list_cluster[order(vgg16_feature_list_cluster[,m]),]
          print("image_size == tile")
        } else {
          print("Error! Define image_size.")
        }
        
        # Separate UMAP coordinates
        rd1 <- vgg16_feature_list_cluster[,c("X1","X2")]
        
        
        # Select labels for plots
        # Modify data
        vgg16_feature_list_cluster <- vgg16_feature_list_cluster %>% dplyr::mutate(
          Diagnosis = ifelse(diagnosis_abbr_sample_HEproject == "PreMDS", "MDS", diagnosis_abbr_sample_HEproject),
          Diagnosis = as.factor(Diagnosis)
        )
        
        
        ################# Convert data into single-cell experiment ##################################################################################################
        
        
        print("Convert data into single-cell experiment")
        # Visual feature matrix
        means <- vgg16_feature_list_cluster[, str_detect(colnames(vgg16_feature_list_cluster), "^V[[:digit:]]*$")] %>% t() %>% as.data.frame()
        
        # Change all values to be ≥0
        ## Although this is needed for some calculations, please use raw values for plotting
        minimal <- sapply(means, function(x) min(x)) %>% min()
        counts <- sapply(means, function(x) x + abs(minimal))
        
        ## Name columns and rows
        rownames(counts) <- paste0('G',1:dim(counts)[1])
        colnames(counts) <- paste0('c',1:dim(counts)[2])
        
        # Make data as single-cell experiment
        sim <- SingleCellExperiment(assays = List(counts = counts))
        
        # Make list of dimension reduction results
        reducedDims(sim)$UMAP <- as.matrix(rd1)
        
        
        ################## Cluster data ##################################################################################################
        
        
        # Define clustering method
        cl1 <- vgg16_feature_list_cluster[,m] # Has to be factor variable
        colData(sim)$KMeans <- cl1
        
        
        # Define dimension reduction and clusters
        sim <- slingshot(sim, clusterLabels = 'KMeans', reducedDim = 'UMAP')
        
        
        ################## Prepare plots with slingshot and colored points ##################################################################################################
        
        
        # # Define colors and make plots
        # 
        # ## If only one pseudotime axis
        # # colors <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[-6])(100)
        # colors <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[1:5])(100)
        # plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
        # 
        # 
        # 
        # ## If 3 pseudotime axes
        # colors0 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[6:4])(100)
        # colors1 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[1:3])(100)
        # colors2 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[8:10])(100)
        # colors3 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Spectral')[6:9])(100)
        # 
        # plotcol0 <- colors0[cut(sim$slingPseudotime_1+sim$slingPseudotime_2+sim$slingPseudotime_3, breaks=100)]
        # plotcol1 <- colors1[cut(sim$slingPseudotime_1, breaks=100)]
        # plotcol2 <- colors2[cut(sim$slingPseudotime_2, breaks=100)]
        # plotcol3 <- colors3[cut(sim$slingPseudotime_3, breaks=100)]
        # 
        # plot(rd1, col = plotcol1, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol2, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol3, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol0, pch=16, asp = 1)
        # lines(SlingshotDataSet(sim), lwd=2, col='black')
        # 
        # 
        # 
        # 
        # ## If 4 pseudotime axes
        # colors0 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Paired')[6:4])(100)
        # colors1 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Paired')[1:3])(100)
        # colors2 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Paired')[7:10])(100)
        # colors3 <- colorRampPalette(brewer.pal(length(summary(cl1)),'Paired')[10:12])(100)
        # 
        # plotcol0 <- colors0[cut(sim$slingPseudotime_1+sim$slingPseudotime_2+sim$slingPseudotime_3+sim$slingPseudotime_4, breaks=100)]
        # plotcol1 <- colors1[cut(sim$slingPseudotime_1, breaks=100)]
        # plotcol2 <- colors2[cut(sim$slingPseudotime_2, breaks=100)]
        # plotcol3 <- colors3[cut(sim$slingPseudotime_3, breaks=100)]
        # plotcol4 <- colors3[cut(sim$slingPseudotime_4, breaks=100)]
        # 
        # # plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
        # plot(rd1, col = plotcol1, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol2, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol3, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol4, pch=16, asp = 1)
        # par(new=TRUE)
        # plot(rd1, col = plotcol0, pch=16, asp = 1)
        # lines(SlingshotDataSet(sim), lwd=2, col='black')
        
        
        ################## Prepare parameters for lineage and slingshot plots ##################################################################################################
        
        
        # We can also see how the lineage structure was intially estimated by the cluster-based minimum spanning tree by using the type argument.
        ## 1. Define lineages
        lin1 <- getLineages(rd1, vgg16_feature_list_cluster[,m])
        
        ## 2. Define lineages + starting  cluster
        lin2 <- getLineages(rd1, vgg16_feature_list_cluster[,m],
                            start.clus = lin1@lineages$Lineage1[1])   # first cluster of first lineage
        
        ## 3. Define both starting and ending clusters
        lin3 <- getLineages(rd1, vgg16_feature_list_cluster[,m],
                            start.clus = lin1@lineages$Lineage1[1],   # first cluster of first lineage
                            end.clus = lin1@lineages$Lineage1[as.numeric(summary(lin1@lineages[which.max(lapply(lin1@lineages, length))])[1])])   # last cluster of longest lineage
        
        # Add slingshot curves
        crv1 <- getCurves(lin1)
        crv1
        
        
        # Legend function
        add_legend <- function(...) {
          opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                      mar=c(0, 0, 0, 0), new=TRUE, mgp=c(2,1,0))
          on.exit(par(opar))
          # plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
          legend(...)
        }
        
        lbl <- c(m, "Diagnosis")
        
        
        ################## Export lineage and slingshot plots ##################################################################################################
        
        
        for (i in lbl) {
          print(i)
          # Plot with dg colors
          if(i == "Diagnosis" & (dg_list == "MDS_CO" | dg_list == "MDS_MDSMPN_AA_CO" | dg_list == "MDS_MDSMPN_CO")) {
            print(paste0("Generate lineage plot 1 for ", i))
            # Add lineages on images
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, ".png"), width = 10, height = 8, units = 'in', res = 300)
            plot(rd1,
                 # bg = brewer.pal(length(summary(i)),"Set1")[i],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size
            points(vgg16_feature_list_cluster[vgg16_feature_list_cluster$Diagnosis=="Healthy",c("X1","X2")],
                   # bg = "red",  # fill color
                   bg = brewer.pal(length(summary(vgg16_feature_list_cluster$Diagnosis)),"Set1")[which(levels(vgg16_feature_list_cluster$Diagnosis)=="Healthy")],
                   col = "black",   # outer color
                   pch = 21,  # point type
                   cex = 2)  # point size
            lines(SlingshotDataSet(sim), lwd=6, type = 'lineages', col = 'black')
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            ## Plot lineage structure
            print(paste0("Generate lineage plot 2 for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, "_startpointgreen.png"), width = 10, height = 8, units = 'in', res = 300)
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[vgg16_feature_list_cluster[,i]],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size)
            points(vgg16_feature_list_cluster[vgg16_feature_list_cluster$Diagnosis=="Healthy",c("X1","X2")],
                   # bg = "red",  # fill color
                   bg = brewer.pal(length(summary(vgg16_feature_list_cluster$Diagnosis)),"Set1")[which(levels(vgg16_feature_list_cluster$Diagnosis)=="Healthy")],
                   col = "black",   # outer color
                   pch = 21,  # point type
                   cex = 2)  # point size
            lines(lin2, lwd = 6, col = 'black', show.constraints = TRUE)
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            ## Plot lineage structure
            print(paste0("Generate lineage plot 3 for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, "_startpointgreen_endpointred.png"), width = 10, height = 8, units = 'in', res = 300)
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[vgg16_feature_list_cluster[,i]],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size)
            points(vgg16_feature_list_cluster[vgg16_feature_list_cluster$Diagnosis=="Healthy",c("X1","X2")],
                   # bg = "red",  # fill color
                   bg = brewer.pal(length(summary(vgg16_feature_list_cluster$Diagnosis)),"Set1")[which(levels(vgg16_feature_list_cluster$Diagnosis)=="Healthy")],
                   col = "black",   # outer color
                   pch = 21,  # point type
                   cex = 2)  # point size
            lines(lin3, lwd = 6, col = 'black', show.constraints = TRUE)
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            print("Generate slingshot plot")
            
            print(paste0("Generate slingshot plot for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Slingshot_", i, "colors.png"), width = 10, height = 8, units = 'in', res = 300)
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[vgg16_feature_list_cluster[,i]],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)
            points(vgg16_feature_list_cluster[vgg16_feature_list_cluster[,i]=="Healthy",c("X1","X2")],
                   # bg = "red",  # fill color
                   bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[which(levels(vgg16_feature_list_cluster[,i])=="Healthy")],
                   col = "black",   # outer color
                   pch = 21,  # point type
                   cex = 2)  # point size
            lines(crv1, lwd = 6, col = 'black')
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16, 
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
          } else {
            print(paste0("Generate lineage plot 1 for ", i))
            # Add lineages on images
            # plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 1)
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, ".png"), width = 10, height = 8, units = 'in', res = 300)
            plot(rd1,
                 # bg = brewer.pal(length(summary(i)),"Set1")[i],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size
            lines(SlingshotDataSet(sim), lwd=6, type = 'lineages', col = 'black')
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            # ## Check lineage structure
            # lin1
            ## Plot lineage structure
            print(paste0("Generate lineage plot 2 for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, "_startpointgreen.png"), width = 10, height = 8, units = 'in', res = 300)
            # par(mgp=c(2,1,0))  # Axis labels closer to axis
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[vgg16_feature_list_cluster[,i]],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size)
            lines(lin2, lwd = 6, col = 'black', show.constraints = TRUE)
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            ## Plot lineage structure
            print(paste0("Generate lineage plot 3 for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Lineageplot_", i, "_startpointgreen_endpointred.png"), width = 10, height = 8, units = 'in', res = 300)
            # par(mgp=c(2,1,0))  # Axis labels closer to axis
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[vgg16_feature_list_cluster[,i]],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)  # point size)
            lines(lin3, lwd = 6, col = 'black', show.constraints = TRUE)
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16,
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
            
            
            
            print(paste0("Generate slingshot plot for ", i))
            png(paste0("PCA_and_UMAP/", y, "/", m, "/Slingshot_", i, "colors.png"), width = 10, height = 8, units = 'in', res = 300)
            # par(mgp=c(2,1,0))  # Axis labels closer to axis
            plot(rd1,
                 # bg = brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1")[cl1],  # fill color
                 bg = c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3"))[vgg16_feature_list_cluster[,i]],  # fill color
                 col = "black",   # outer color
                 # asp = 1,  # dimension control
                 pch = 21,  # point type
                 cex = 2,  # point size
                 cex.lab = 0.85,
                 xlim = c(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))+2.5),
                 ylim = c(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))),
                 axes = FALSE)
            lines(crv1, lwd = 6, col = 'black')
            add_legend(x = max(ceiling(vgg16_feature_list_cluster$X1)), y= max(vgg16_feature_list_cluster$X2),
                   legend=unique(vgg16_feature_list_cluster[,i]),
                   pch=16, 
                   col = unique(c(brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set1"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set2"), brewer.pal(length(summary(vgg16_feature_list_cluster[,i])),"Set3")))[unique(vgg16_feature_list_cluster[,i])],  # fill color
                   # col = Set4,  # fill color
                   horiz=FALSE,
                   cex=1)
            axis(1,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X1)), seq(min(floor(vgg16_feature_list_cluster$X1)), max(ceiling(vgg16_feature_list_cluster$X1))), max(ceiling(vgg16_feature_list_cluster$X1))))
            axis(2,cex.axis=0.75, at = c(min(floor(vgg16_feature_list_cluster$X2)), seq(min(floor(vgg16_feature_list_cluster$X2)), max(ceiling(vgg16_feature_list_cluster$X2))), max(ceiling(vgg16_feature_list_cluster$X2))))
            dev.off()
          }
        }
        
      }
    }
  }
}
