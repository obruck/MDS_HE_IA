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
input_directory = "DEFINE WHERE YOUR H&E IMAGES ARE LOCATED"
output_directory = 'DEFINE WHERE YOU WANT THE MODIFIED H&E IMAGES TO BE EXPORTED'
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














# Module help
## Which modules are available in library cv2
help(cv2)

## Inspect one module
import types
print([getattr(cv2, a) for a in dir(cv2)
  if isinstance(getattr(cv2, a), types.FunctionType)])

## Inspect one module    
from inspect import getmembers, isfunction
functions_list = [o for o in getmembers(cv2) if isfunction(o[1])]