#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:34:28 2019

@author: oscarbruck

# Convert RGB images grayscale
"""


# Import modules
import os
import cv2


# Define parameters
input_directory = "/data/images/1_preprocessed_tma_color/"
output_directory = "/data/images/2_preprocessed_tma_grayscale/"


# Convert colors
for filename in os.listdir(input_directory):
    if filename.endswith(".jpg"):
        #print(os.path.join(input_directory, filename))
        img_original = cv2.imread(input_directory+filename, 1)
        img = cv2.cvtColor(img_original, cv2.COLOR_BGR2GRAY)        #make as grayscale
        output_file_name = output_directory+filename
        cv2.imwrite(output_file_name, img)
    else:
        break

cv2.destroyAllWindows()
