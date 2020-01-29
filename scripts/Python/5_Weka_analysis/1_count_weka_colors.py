#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:17:28 2019

@author: oscarbruck

# Script to detect and dissect colors and calculate their pixel amount in regards to all other pixels
"""




# Import modules
import numpy as np
import cv2
import os
import pandas as pd


# Define parameters
input_directory = "/data/images/4_weka_images/"
output_directory = "/data/processed_files/"

# Work out what we are looking for
## Bone marrow stainings are here dissected into four classes.
## Colors have been defined based on values exported from Weka analysis.
background_pixel = [198,118,255]
stroma_pixel = [79,255,130]
tumor_pixel = [255,0,0]
rbc_pixel = [255,225,10]

# Define variables you need to export
background = list()
stroma = list()
tumor = list()
rbc = list()
imagelist = list()
totalpixel_wo_bg = list()
background_prop = list()
result_totalpixel_wo_bg = []
stroma_prop = list()
tumor_prop = list()
rbc_prop = list()
df = []
final = pd.DataFrame()




# Looping function
for filename in os.listdir(input_directory):
    if filename.endswith(".png"):
        image = cv2.imread(input_directory+filename)
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        # Resize image = faster
        modified_image = cv2.resize(image, (600, 400), interpolation = cv2.INTER_AREA)
        modified_image = modified_image.reshape(modified_image.shape[0]*modified_image.shape[1], 3)
        # Find all pixels where the 3 RGB values match "background_pixel" etc, and count
        result_bg = np.count_nonzero(np.all(modified_image==background_pixel,axis=1))
        result_st = np.count_nonzero(np.all(modified_image==stroma_pixel,axis=1))
        result_tu = np.count_nonzero(np.all(modified_image==tumor_pixel,axis=1))
        result_rbc = np.count_nonzero(np.all(modified_image==rbc_pixel,axis=1))
        # Append results
        background.append(result_bg)
        stroma.append(result_st)
        tumor.append(result_tu)
        rbc.append(result_rbc)
        imagelist.append(filename)
        # Calculate proportions
        result_totalpixel = result_bg+result_st+result_tu+result_rbc
        result_totalpixel_wo_bg = result_st+result_tu+result_rbc
        totalpixel_wo_bg.append(result_totalpixel_wo_bg)
        result_background_prop = result_bg/result_totalpixel
        background_prop.append(result_background_prop)
        if result_totalpixel_wo_bg > 0:
            result_stroma_prop = result_st/result_totalpixel_wo_bg
            result_tumor_prop = result_tu/result_totalpixel_wo_bg
            result_rbc_prop = result_rbc/result_totalpixel_wo_bg
            # Combine results into a list
            # df = np.column_stack((background,imagelist))
            # Combine results into a dataframe
        else:
            result_stroma_prop = 0
            result_tumor_prop = 0
            result_rbc_prop = 0
        rbc_prop.append(result_rbc_prop)
        stroma_prop.append(result_stroma_prop)
        tumor_prop.append(result_tumor_prop)
        print(len(tumor_prop))
        final = pd.DataFrame(
                {'background': background,
                 'stroma': stroma,
                 'tumor': tumor,
                 'rbc': rbc,
                 'totalpixel_wo_bg': totalpixel_wo_bg,
                 'background_proportion': background_prop,
                 'stroma_proportion': stroma_prop,
                 'tumor_proportion': tumor_prop,
                 'rbc_proportion': rbc_prop,
                 'filename': imagelist,
                 })
    else:
        pass

final.to_excel(output_directory+'Color_proportions_tilelevel.xlsx', sheet_name='sheet1', index=False)


