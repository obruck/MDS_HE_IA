#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 08:43:19 2019

@author: oscarbruck

# Uniform background color around tissues
"""

# Import modules
import os
import cv2
import numpy as np
from skimage import util

# Define parameters
input_directory = "/data/images/0_raw/"
output_directory = '/data/images/1_preprocessed_tma_color/'

# Remove background
for filename in os.listdir(input_directory):
    if filename.endswith(".jpeg"):
        img_original = cv2.imread(input_directory+filename)
        img = cv2.cvtColor(img_original, cv2.COLOR_BGR2GRAY)
        thresh = 200
        img = cv2.threshold(img, thresh, 255, cv2.THRESH_BINARY)[1]
        img = util.invert(img)
        #make as binary
        img[img!=0] = 255
        # dilate images
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (20, 20))
        img = cv2.dilate(img, kernel)
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
        img_final = cv2.bitwise_and(img_original, img)
        img_final[np.where((img_final==[0,0,0]).all(axis=2))] = [255,255,255]
        output_file_name = output_directory+filename+"_modified.jpg"
        cv2.imwrite(output_file_name, img_final)
    else:
        break


cv2.destroyAllWindows()
