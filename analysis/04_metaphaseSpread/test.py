from skimage.filters import threshold_local, sobel
from skimage.morphology import binary_dilation, binary_erosion, dilation
from skimage.measure import label, regionprops
import shared.objects as obj
from skimage.segmentation import watershed
import shared.segmentation as seg
import skimage.io as skio
from skimage.morphology import extrema
import shared.image as ima
from scipy import ndimage
import napari
import shared.display as dis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220909_EdU_metaphasetest/20220909_EdU_metaphasetest/"
start_fov = 21
total_fov = 26

for f in range(total_fov):
    fov = f + start_fov
    file_name = 'DM_%s_RAW' % fov
    img_hoechst = skio.imread("%s%s_ch00.tif" % (master_folder, file_name), plugin="tifffile")
    img_FISH = skio.imread("%s%s_ch01.tif" % (master_folder, file_name), plugin="tifffile")
    img_EdU = skio.imread("%s%s_ch02.tif" % (master_folder, file_name), plugin="tifffile")

    viewer1 = napari.Viewer()
    viewer1.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer1.add_image(img_FISH, blending='additive', colormap='green', contrast_limits=[0, img_FISH.max()])
    viewer1.add_image(img_EdU, blending='additive', colormap='magenta', contrast_limits=[0, 10000])
    napari.run()
