# Created: 31.03.2019, OB

# Build elastic net models to predict TET2 mutation using feature matrix extracted using

## VGG16 network --> lines 9-142
## Xception network --> lines 143-273


############################### VGG16 - Load packages ################################

library(tidyverse)
library(dplyr)
library(glmnet)
library(readr)
library(stringr)
library(lubridate)
library(reshape2)
library(parallel)
library(readxl)
library(biglasso)
library(doParallel)

# Create cluster
nodes <- detectCores()-1
print(nodes)
cl <- makeCluster(nodes)
doParallel::registerDoParallel(nodes)


############################### VGG16 - Load data ##############################################################


# Load feature matrix data
vgg16_feature_list1 <- readRDS("/data/processed_files/Feature_list_VGG16_MDS.rds")

## Make a variable with the image spot name
vgg16_feature_list1$name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles/", "", vgg16_feature_list1$vgg16_feature_list_name)
vgg16_feature_list1$name <- gsub("\\/MDS_TMA", "MDS_TMA", vgg16_feature_list1$name)
vgg16_feature_list1$name <- gsub("_modified[[:print:]]+", "", vgg16_feature_list1$name)
vgg16_feature_list1 <- vgg16_feature_list1 %>% dplyr::select(name, everything())

# Load clinical data
vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")

# Keep only MDS patients
vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample=="MDS" | vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample=="PreMDS",]
## Modify data
hm1 <- vgg16_feature_list1 %>% filter(name %in% vgg16_feature_list_cluster_clinical$Image_name)

## Merge
clinical <- vgg16_feature_list_cluster_clinical %>%
  dplyr::select(Image_name, TET2_kval) %>%
  rename(name = Image_name, 
         cmplx_ipss = TET2_kval)
hm1 <- clinical %>%
  left_join(vgg16_feature_list1, by = "name") %>%
  dplyr::select(name, vgg16_feature_list_name, cmplx_ipss, everything())


########################### VGG16 - Create train dataset ############################################################


#Glmnet
set.seed(101)
## Dg classes
hm1$cmplx_ipss <- as.numeric(hm1$cmplx_ipss)

## Remove NA
hm1 <- hm1[!is.na(hm1$cmplx_ipss),]
hm1 <- hm1[!is.na(hm1$vgg16_feature_list_name),]


## Sample data to training and test set
### Sample by wsi number
n <- length(unique(hm1$name))
sample <- sample(unique(hm1$name), size = n * 0.67, replace = FALSE)


### TrainX
train_mds <- hm1 %>% filter(name %in% sample)
trainX_mds <- train_mds[,3:ncol(hm1)] %>% dplyr::select(-cmplx_ipss)
trainX_mds <- sapply(trainX_mds, function(x) as.numeric(x))
trainX_mds <- as.matrix(trainX_mds)


### TrainY
trainY_mds <- as.factor(as.matrix(train_mds["cmplx_ipss"]))

# Make biglasso matric
trainX_mds.bm <- as.big.matrix(trainX_mds)


############################# VGG16 - Elastic net with biglasso #############################################


# Big lasso for elastic net
set.seed(101)

# Big lasso CV
fit.bin.lasso_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "lasso", alpha = 1)
fit.bin.ridge_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "ridge", alpha = 0)
fit.bin.enet_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "enet", alpha = 0.5)

# Lambda.1se
fit.bin.lasso_min_cv$cve_1se <- fit.bin.lasso_min_cv$cve[fit.bin.lasso_min_cv$lambda == fit.bin.lasso_min_cv$lambda.min] + fit.bin.lasso_min_cv$cvse[fit.bin.lasso_min_cv$lambda == fit.bin.lasso_min_cv$lambda.min]
fit.bin.lasso_min_cv$lambda.1se <- fit.bin.lasso_min_cv$lambda[which.min(abs(fit.bin.lasso_min_cv$cve_1se-fit.bin.lasso_min_cv$cve))]

fit.bin.ridge_min_cv$cve_1se <- fit.bin.ridge_min_cv$cve[fit.bin.ridge_min_cv$lambda == fit.bin.ridge_min_cv$lambda.min] + fit.bin.ridge_min_cv$cvse[fit.bin.ridge_min_cv$lambda == fit.bin.ridge_min_cv$lambda.min]
fit.bin.ridge_min_cv$lambda.1se <- fit.bin.ridge_min_cv$lambda[which.min(abs(fit.bin.ridge_min_cv$cve_1se-fit.bin.ridge_min_cv$cve))]

fit.bin.enet_min_cv$cve_1se <- fit.bin.enet_min_cv$cve[fit.bin.enet_min_cv$lambda == fit.bin.enet_min_cv$lambda.min] + fit.bin.enet_min_cv$cvse[fit.bin.enet_min_cv$lambda == fit.bin.enet_min_cv$lambda.min]
fit.bin.enet_min_cv$lambda.1se <- fit.bin.enet_min_cv$lambda[which.min(abs(fit.bin.enet_min_cv$cve_1se-fit.bin.enet_min_cv$cve))]


# Fit values according to lambda.min & lambda.1se and alpha = 0, 1 and enet
fit.bin.lasso_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='lasso', alpha = 1, family = "binomial", lambda = fit.bin.lasso_min_cv$lambda.min)
fit.bin.ridge_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='ridge', alpha = 0, family = "binomial", lambda = fit.bin.ridge_min_cv$lambda.min)
fit.bin.enet_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='enet', alpha = 0.5, family = "binomial", lambda = fit.bin.enet_min_cv$lambda.min)
fit.bin.lasso_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='lasso', alpha = 1, family = "binomial", lambda = fit.bin.lasso_min_cv$lambda.1se)
fit.bin.ridge_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='ridge', alpha = 0, family = "binomial", lambda = fit.bin.ridge_min_cv$lambda.1se)
fit.bin.enet_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='enet', alpha = 0.5, family = "binomial", lambda = fit.bin.enet_min_cv$lambda.1se)

# Generate predictions
pred_lasso_min <- predict(fit.bin.lasso_min, trainX_mds.bm, type="response")
pred_ridge_min <- predict(fit.bin.ridge_min, trainX_mds.bm, type="response")
pred_enet_min <- predict(fit.bin.enet_min, trainX_mds.bm, type="response")
pred_lasso_1se <- predict(fit.bin.lasso_1se, trainX_mds.bm, type="response")
pred_ridge_1se <- predict(fit.bin.ridge_1se, trainX_mds.bm, type="response")
pred_enet_1se <- predict(fit.bin.enet_1se, trainX_mds.bm, type="response")



# Save list of glmnet models and optimla cv parameters (alpha and lambda)
#cv_mds_list <- list(cv1, md_mds_lambda1se, md_mds_lambdamin, fit_mds_l1, fit_mds_l2, cv_mds_l1_lambda1se, cv_mds_l2_lambda1se)
cv_mds_list <- list(fit.bin.enet_min_cv, fit.bin.enet_min, fit.bin.enet_1se, fit.bin.lasso_min_cv, fit.bin.lasso_min, fit.bin.lasso_1se, fit.bin.ridge_min_cv, fit.bin.ridge_min, fit.bin.ridge_1se, pred_lasso_min, pred_ridge_min, pred_enet_min, pred_lasso_1se, pred_ridge_1se, pred_enet_1se) 
saveRDS(cv_mds_list, "/data/processed_files/Elastic_net_results/tet2/VGG16/vgg16_pred_tet2_elastinet_biglasso_MDSonly.rds")

stopCluster(cl)

rm(list=ls())


############################### Xception - Load packages ##############################################################


library(tidyverse)
library(dplyr)
library(glmnet)
library(readr)
library(stringr)
library(lubridate)
library(reshape2)
library(parallel)
library(readxl)
library(biglasso)
library(doParallel)

# Create cluster
nodes <- detectCores()-1
#nodes <- 20
print(nodes)
cl <- makeCluster(nodes)
doParallel::registerDoParallel(nodes)


########################### Xception - Load data ##############################################################


# Load feature matrix data
vgg16_feature_list1 <- readRDS("/data/processed_files/Feature_list_VGG16_MDS.rds")

## Make a variable with the image spot name
vgg16_feature_list1$name <- gsub("/csc/mustjoki/imageanalysis/analysisR/input/mds_he_tiles/", "", vgg16_feature_list1$vgg16_feature_list_name)
vgg16_feature_list1$name <- gsub("\\/MDS_TMA", "MDS_TMA", vgg16_feature_list1$name)
vgg16_feature_list1$name <- gsub("_modified[[:print:]]+", "", vgg16_feature_list1$name)
vgg16_feature_list1 <- vgg16_feature_list1 %>% dplyr::select(name, everything())

# Load clinical data
vgg16_feature_list_cluster_clinical <- read_xlsx("/data/tables/Clinicaldata_forimages.xlsx")

# Keep only MDS patients
vgg16_feature_list_cluster_clinical <- vgg16_feature_list_cluster_clinical[vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample=="MDS" | vgg16_feature_list_cluster_clinical$diagnosis_abbr_sample=="PreMDS",]
## Modify data
hm1 <- vgg16_feature_list1 %>% filter(name %in% vgg16_feature_list_cluster_clinical$Image_name)

## Merge
clinical <- vgg16_feature_list_cluster_clinical %>%
  dplyr::select(Image_name, TET2_kval) %>%
  rename(name = Image_name,
         cmplx_ipss = TET2_kval)
hm1 <- clinical %>%
  left_join(vgg16_feature_list1, by = "name") %>%
  dplyr::select(name, vgg16_feature_list_name, cmplx_ipss, everything())


########################### Xception - Create train dataset ############################################################


#Glmnet
set.seed(101)
## Dg classes
hm1$cmplx_ipss <- as.numeric(hm1$cmplx_ipss)

## Remove NA
hm1 <- hm1[!is.na(hm1$cmplx_ipss),]
hm1 <- hm1[!is.na(hm1$vgg16_feature_list_name),]


## Sample data to training and test set
### Sample by wsi number
n <- length(unique(hm1$name))
sample <- sample(unique(hm1$name), size = n * 0.67, replace = FALSE)


### TrainX
train_mds <- hm1 %>% filter(name %in% sample)
trainX_mds <- train_mds[,3:ncol(hm1)] %>% dplyr::select(-cmplx_ipss)
trainX_mds <- sapply(trainX_mds, function(x) as.numeric(x))
trainX_mds <- as.matrix(trainX_mds)

### TrainY
trainY_mds <- as.factor(as.matrix(train_mds["cmplx_ipss"]))

# Make biglasso matric
trainX_mds.bm <- as.big.matrix(trainX_mds)


########################### Xception - Elastic net with biglasso ############################################################


# Big lasso for elastic net
set.seed(101)

# Big lasso CV
fit.bin.lasso_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "lasso", alpha = 1)
fit.bin.ridge_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "ridge", alpha = 0)
fit.bin.enet_min_cv <- cv.biglasso(trainX_mds.bm, trainY_mds, family = "binomial", nfolds = 5, seed = 123, penalty = "enet", alpha = 0.5)

# Lambda.1se
fit.bin.lasso_min_cv$cve_1se <- fit.bin.lasso_min_cv$cve[fit.bin.lasso_min_cv$lambda == fit.bin.lasso_min_cv$lambda.min] + fit.bin.lasso_min_cv$cvse[fit.bin.lasso_min_cv$lambda == fit.bin.lasso_min_cv$lambda.min]
fit.bin.lasso_min_cv$lambda.1se <- fit.bin.lasso_min_cv$lambda[which.min(abs(fit.bin.lasso_min_cv$cve_1se-fit.bin.lasso_min_cv$cve))]

fit.bin.ridge_min_cv$cve_1se <- fit.bin.ridge_min_cv$cve[fit.bin.ridge_min_cv$lambda == fit.bin.ridge_min_cv$lambda.min] + fit.bin.ridge_min_cv$cvse[fit.bin.ridge_min_cv$lambda == fit.bin.ridge_min_cv$lambda.min]
fit.bin.ridge_min_cv$lambda.1se <- fit.bin.ridge_min_cv$lambda[which.min(abs(fit.bin.ridge_min_cv$cve_1se-fit.bin.ridge_min_cv$cve))]

fit.bin.enet_min_cv$cve_1se <- fit.bin.enet_min_cv$cve[fit.bin.enet_min_cv$lambda == fit.bin.enet_min_cv$lambda.min] + fit.bin.enet_min_cv$cvse[fit.bin.enet_min_cv$lambda == fit.bin.enet_min_cv$lambda.min]
fit.bin.enet_min_cv$lambda.1se <- fit.bin.enet_min_cv$lambda[which.min(abs(fit.bin.enet_min_cv$cve_1se-fit.bin.enet_min_cv$cve))]


# Fit values according to lambda.min & lambda.1se and alpha = 0, 1 and enet
fit.bin.lasso_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='lasso', alpha = 1, family = "binomial", lambda = fit.bin.lasso_min_cv$lambda.min)
fit.bin.ridge_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='ridge', alpha = 0, family = "binomial", lambda = fit.bin.ridge_min_cv$lambda.min)
fit.bin.enet_min <- biglasso(trainX_mds.bm, trainY_mds, penalty ='enet', alpha = 0.5, family = "binomial", lambda = fit.bin.enet_min_cv$lambda.min)
fit.bin.lasso_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='lasso', alpha = 1, family = "binomial", lambda = fit.bin.lasso_min_cv$lambda.1se)
fit.bin.ridge_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='ridge', alpha = 0, family = "binomial", lambda = fit.bin.ridge_min_cv$lambda.1se)
fit.bin.enet_1se <- biglasso(trainX_mds.bm, trainY_mds, penalty ='enet', alpha = 0.5, family = "binomial", lambda = fit.bin.enet_min_cv$lambda.1se)

# Generate predictions
pred_lasso_min <- predict(fit.bin.lasso_min, trainX_mds.bm, type="response")
pred_ridge_min <- predict(fit.bin.ridge_min, trainX_mds.bm, type="response")
pred_enet_min <- predict(fit.bin.enet_min, trainX_mds.bm, type="response")
pred_lasso_1se <- predict(fit.bin.lasso_1se, trainX_mds.bm, type="response")
pred_ridge_1se <- predict(fit.bin.ridge_1se, trainX_mds.bm, type="response")
pred_enet_1se <- predict(fit.bin.enet_1se, trainX_mds.bm, type="response")


# Save list of glmnet models and optimla cv parameters (alpha and lambda)
#cv_mds_list <- list(cv1, md_mds_lambda1se, md_mds_lambdamin, fit_mds_l1, fit_mds_l2, cv_mds_l1_lambda1se, cv_mds_l2_lambda1se)
cv_mds_list <- list(fit.bin.enet_min_cv, fit.bin.enet_min, fit.bin.enet_1se, fit.bin.lasso_min_cv, fit.bin.lasso_min, fit.bin.lasso_1se, fit.bin.ridge_min_cv, fit.bin.ridge_min, fit.bin.ridge_1se, pred_lasso_min, pred_ridge_min, pred_enet_min, pred_lasso_1se, pred_ridge_1se, pred_enet_1se) 
saveRDS(cv_mds_list, "/data/processed_files/Elastic_net_results/tet2/Xception/xception_pred_tet2_elastinet_biglasso_MDSonly.rds")

stopCluster(cl)
