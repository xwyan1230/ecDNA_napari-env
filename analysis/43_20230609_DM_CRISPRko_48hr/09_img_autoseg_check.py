import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import math
import shared.objects as obj
import tifffile as tif
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'C'
sample = 'C3'
# batch = 2
total_fov = 16
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
convex_conversion_threshold = 0.85
local_size = 200
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi

for fov in range(total_fov):
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_4pos_s%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 11000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 20)

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    napari.run()

print("DONE!")