import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'dWee1_125nM_24hr'
file_name = '20240307_sp8_iWee1anddWee1_acidFISH_24hr_dWee1_125nM_24hr_9pos_1'
total_fov = 9

for i in range(total_fov):
    img_hoechst = skio.imread("%s%s/%s_s%s_ch00.tif" % (data_dir, sample, file_name, i), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_s%s_ch01.tif" % (data_dir, sample, file_name, i), plugin="tifffile")
    # img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    # img_DNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    napari.run()