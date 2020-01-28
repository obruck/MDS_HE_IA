# Created: 03.01.2020, OB

# Script to extract visual features with the Xception network


######################## LOAD PACKAGES ##################################################################


library(reticulate)
use_virtualenv("/home/jimm/anaconda3/envs/rpyenv")
library(keras)
library(tidyverse)
library(readxl)
library(lubridate)
library(data.table)


######################## MODEL ##################################################################


model <- keras::application_xception(weights = "imagenet", include_top = FALSE, pooling = "avg")
model


###################### FUNCTION ######################################################################


# Load data

# Helper function for reading in images and preprocessing them
image_prep <- function(x) {
arrays <- lapply(x, function(path) {
img <- image_load(path, target_size = c(299, 299))
x <- image_to_array(img)
x <- array_reshape(x, c(1, dim(x)))
x <- xception_preprocess_input(x)
})
do.call(abind::abind, c(arrays, list(along = 1)))
}


###################### EXTRACT FEATURES ######################################################################


# Directory path
## Don't use "/" at the end of the path
image_files_path <- "WHERE_YOUR_PREPROCESSED_IMAGES_ARE_LOCATED"
output_directory <- "WHERE_YOUR_IMAGE_FEATURES_WILL_BE_EXPORTED"
dir_list <- list.dirs(image_files_path, full.names = TRUE)
dir_list <- base::sample(list.dirs(image_files_path, full.names = TRUE), 5)

# Extract features by looping each image of each directory
## looping group = each directory
for (image_dir in dir_list) {
  print(paste0("Extracting images from ", image_dir))
  image_dir_name <- gsub("[[:print:]]*\\/", "", image_dir)
  print(paste0("Saving as ", image_dir_name))
  image_list <- list.files(image_dir, full.names = TRUE, recursive = TRUE)
  image_list <- purrr::discard(image_list,.p = ~stringr::str_detect(.x,"TXT"))
  vgg16_feature_list <- data.frame()
  for (image in image_list) {
    print(image)
    cat("Image", which(image_list == image), "from", length(image_list))
    
    vgg16_feature <- predict(model, image_prep(image))
    
    flatten <- as.data.frame.table(vgg16_feature, responseName = "value") %>%
      dplyr::select(value)
    flatten <- cbind(image, as.data.frame(t(flatten)))
    
    vgg16_feature_list <- rbind(vgg16_feature_list, flatten)
    
  }
  save(vgg16_feature_list, file = paste0(output_directory, "/MDS_imagetiles_xception_", image_dir_name, ".Rdata"))
  rm(vgg16_feature_list)
}


###################### BIND TILES INTO ONE LARGE FEATURE MATRIX ######################################################################


# Load data
## Define filenames
setwd(output_directory)
files <- list.files(pattern = ".Rdata")
## Read and append data files as a list
datalist = list()
for (i in 1:length(files)) {
  load(files[i])
  print(paste0("Loading ", i))
  datalist[[i]] <- vgg16_feature_list
}

## Merge datafiles together
feature_list = do.call(rbind, datalist)


###################### SAVE DATA ######################################################################


# Save feature list
saveRDS(feature_list, "Feature_list_xception_MDS.rds")
