from skimage.measure import regionprops, label
import shared.image as ima
import skimage.io as skio
import numpy as np
import shared.dataframe as dat
import pandas as pd
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220927_Jun_EGFR_RPAs33p_Edu/EGFR_RPAs33p_Edu/"
sample = 'gbm39ec con'
master_path = '%s%s/' % (master_folder, sample)
end_fov = 10
start_fov = 1
total_fov = end_fov - start_fov + 1
local_size = 150

# IMAGING ANALYSIS
for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, start nuclear segmentation FOV %s/%s" % (sample, fov, total_fov))
    file_prefix = "40x %s-%s" % (sample, fov)
    # LOAD IMAGE
    img_nuclear = skio.imread("%sC2-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_DNAFISH = skio.imread("%sC1-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_RPA = skio.imread("%sC3-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_EdU = skio.imread("%sC4-%s.tif" % (master_path, file_prefix), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
    # viewer.add_image(img_EdU, blending='additive', colormap='magenta', contrast_limits=[0, img_EdU.max()])
    napari.run()