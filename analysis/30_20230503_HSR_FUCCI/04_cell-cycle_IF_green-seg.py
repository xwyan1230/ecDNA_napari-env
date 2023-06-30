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
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230517_analysis_HSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'HSR_81pos_merge'

# set parameters
pixel_size = 142  # nm (sp8 confocal 40x 2048*2048)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.8, 1.5]  # used to filter nucleus, DM [0.8, 1.3], HSR [0.8, 1.5]
convex_conversion_threshold = 0.85
circ_threshold = 0.8
otsu_factor = 1.0

# segmentation
local_factor_nuclear = 77  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

img_hoechst = skio.imread("%s%s/%s_RAW_ch03.tif" % (data_dir, sample, sample), plugin="tifffile")
img_red = skio.imread("%s%s/%s_RAW_ch00.tif" % (data_dir, sample, sample), plugin="tifffile")
img_green = skio.imread("%s%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")
img_MYC = skio.imread("%s%s/%s_RAW_ch02.tif" % (data_dir, sample, sample), plugin="tifffile")

print("running green seg...")
img_green_seg = seg.nuclear_seg2_print(img_green, local_factor=local_factor_nuclear, bg_thresh=10, min_size=min_size_nuclear, max_size=2*max_size_nuclear)
print("running green seg convex...")
img_green_seg_convex = seg.obj_to_convex_filter_print(img_green_seg, threshold=convex_conversion_threshold)
print("saving green seg convex...")
tif.imwrite("%s%s/HSR_IF_green_seg_convex.tif" % (output_dir, sample), img_green_seg_convex)

print("DONE!")


