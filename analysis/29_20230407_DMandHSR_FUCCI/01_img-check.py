import skimage.io as skio
import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'

img_hoechst = skio.imread("%s%s/%s_RAW_ch03.tif" % (data_dir, sample, sample), plugin="tifffile")[:22000, :33601]
img_red = skio.imread("%s%s/%s_RAW_ch00.tif" % (data_dir, sample, sample), plugin="tifffile")[:22000, :33601]
img_green = skio.imread("%s%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")[:22000, :33601]
img_MYC = skio.imread("%s%s/%s_RAW_ch02.tif" % (data_dir, sample, sample), plugin="tifffile")[:22000, :33601]
# tif.imwrite("%s%s/DM_hoechst_IF_22000-33601.tif" % (output_dir, sample), img_hoechst)
tif.imwrite("%s%s/DM_red_IF_22000-33601.tif" % (output_dir, sample), img_red)
tif.imwrite("%s%s/DM_green_IF_22000-33601.tif" % (output_dir, sample), img_green)
tif.imwrite("%s%s/DM_MYC_IF_22000-33601.tif" % (output_dir, sample), img_MYC)

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
# plt.imsave("%s%s/DM_IF_22000-33601.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()


