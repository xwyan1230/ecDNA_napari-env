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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'

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

img_hoechst = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_red = skio.imread("%s%s/DM_red_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_green = skio.imread("%s%s/DM_green_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_MYC = skio.imread("%s%s/DM_MYC_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")

print("running nuclear seg...")
img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear, max_size=max_size_nuclear)
print("saving nuclear seg...")
tif.imwrite("%s%s/DM_hoechst_IF_22000-33601_seg.tif" % (output_dir, sample), img_nuclear_seg)
img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg, threshold=convex_conversion_threshold))
print("saving nuclear seg...")
tif.imwrite("%s%s/DM_hoechst_IF_22000-33601_seg_convex.tif" % (output_dir, sample), img_nuclear_seg_convex)
print("running red seg...")
img_red_seg = seg.nuclear_seg2(img_red, local_factor=local_factor_nuclear, bg_thresh=10, min_size=min_size_nuclear, max_size=2*max_size_nuclear)
img_red_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_red_seg, threshold=convex_conversion_threshold))
print("saving red seg...")
tif.imwrite("%s%s/DM_red_IF_22000-33601_seg.tif" % (output_dir, sample), img_red_seg_convex)
print("running green seg...")
img_green_seg = seg.nuclear_seg2(img_green, local_factor=local_factor_nuclear, bg_thresh=10, min_size=min_size_nuclear, max_size=2*max_size_nuclear)
img_green_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_green_seg, threshold=convex_conversion_threshold))
print("saving green seg...")
tif.imwrite("%s%s/DM_green_IF_22000-33601_seg.tif" % (output_dir, sample), img_green_seg_convex)

"""viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
viewer.add_image(img_nuclear_seg_convex, blending='additive')
viewer.add_image(img_red_seg_convex, blending='additive')
viewer.add_image(img_green_seg_convex, blending='additive')
# plt.imsave("%s%s/DM_IF_22000-33601.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()"""


