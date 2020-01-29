#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Apr  5 09:07:49 2019

@author: oscarbruck

# Crop TMA images into 500 equally-sized tiles and discard all-white (empty) tiles
"""

# Import modules
import image_slicer
import os
import cv2
import imageio




# Define parameters
input_directory = "/data/images/2_preprocessed_tma_grayscale/"
output_directory = '/data/images/3_preprocessed_tiles_grayscale/'
num_tiles = 500

# Slice images
for filename in os.listdir(input_directory):
        if filename.endswith(".jpg"): 
            tiles = image_slicer.slice(input_directory+filename, num_tiles, save=False)
            image_slicer.save_tiles(tiles, format='jpeg', prefix=filename, directory=output_directory)
            continue
        else:
            pass

# Discard all-white images
for filename in os.listdir(output_directory):
    if filename.endswith(".jpg"): 
        img = imageio.imread(output_directory+filename)
        #img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        #print(filename)
        average = img.mean(axis=0).mean(axis=0)
        if average > 240:
            os.remove(output_directory+filename)
            continue
        else:
            pass
        continue
    else:
        pass

