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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'

# set parameters
pixel_size = 767  # nm (keyence 10x)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.8, 1.5]  # used to filter nucleus, DM [0.8, 1.3], HSR [0.8, 1.5]
convex_conversion_threshold = 0.85
circ_threshold = 0.8
otsu_factor = 1.0

# segmentation
local_factor_nuclear = 77  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
print(min_size_nuclear)
print(max_size_nuclear)

img_hoechst_Keyence = skio.imread("%s%s/Hoechst_cut.tif" % (data_dir, sample), plugin="tifffile")

print("running nuclear seg...")
img_nuclear_seg = seg.nuclear_seg3(img_hoechst_Keyence, min_size=min_size_nuclear, max_size=max_size_nuclear, maxima_threshold=1)
img_nuclear_seg_convex = seg.obj_to_convex_filter_print(img_nuclear_seg, threshold=convex_conversion_threshold)
if not os.path.exists("%s%s/seg_keyence/" % (output_dir, sample)):
    os.makedirs("%s%s/seg_keyence/" % (output_dir, sample))
tif.imwrite("%s%s/seg_keyence/Hoechst_keyence_seg_convex.tif" % (output_dir, sample), img_nuclear_seg_convex)

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_nuclear_seg_convex, blending='additive')
# plt.imsave("%s%s/DM_IF_22000-33601.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()


