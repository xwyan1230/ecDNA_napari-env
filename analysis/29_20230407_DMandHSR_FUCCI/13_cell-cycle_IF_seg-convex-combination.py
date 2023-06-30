import skimage.io as skio
import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import math
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
from skimage.measure import label, regionprops_table, regionprops
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'

img_hoechst = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_nuclear_seg = skio.imread("%s%s/DM_hoechst_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")
img_red_seg = skio.imread("%s%s/DM_red_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")
img_green_seg = skio.imread("%s%s/DM_green_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")

nuclear_red_combination = seg.seg_combination(img_nuclear_seg, img_red_seg)
tif.imwrite("%s%s/DM_nuclear_red_IF_22000-33601_seg_convex.tif" % (output_dir, sample), nuclear_red_combination)
nuclear_red_green_combination = seg.seg_combination(nuclear_red_combination, img_green_seg)
tif.imwrite("%s%s/DM_nuclear_red_green_IF_22000-33601_seg_convex.tif" % (output_dir, sample), nuclear_red_green_combination)


"""viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive')
viewer.add_image(img_nuclear_seg, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_red_seg, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_green_seg, blending='additive', colormap='green', contrast_limits=[0, 65535])
# plt.imsave("%s%s/DM_IF_22000-33601.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()"""


