import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import tifffile as tif
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

sample = 'DM-Ctrl_mix_mCh-BRD4'
batch = 0
if batch > 0:
    data_dir = '%s%s/grid/%s/' % (data_dir, sample, batch)
else:
    data_dir = '%s%s/grid/' % (data_dir, sample)

output_dir = data_dir

start_fov = 0
end_fov = 0
total_fov = end_fov-start_fov+1
group = 'neg'

for f in range(total_fov):
    fov = f+start_fov
    file_name = '%s_%s' % (sample, group)
    img_hoechst = skio.imread("%s%s_hoechst_%s.tif" % (data_dir, file_name, fov), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s_DNAFISH_%s.tif" % (data_dir, file_name, fov), plugin="tifffile")
    img_seg = skio.imread("%s%s_seg_%s.tif" % (data_dir, file_name, fov), plugin="tifffile")
    img_ecseg = skio.imread("%s%s_ecseg1_%s.tif" % (data_dir, file_name, fov), plugin="tifffile")
    img_DNAFISH_filtered = img_DNAFISH.copy()
    img_DNAFISH_filtered[img_ecseg == 1] = 0

    tif.imwrite("%s%s_DNAFISH-filtered1_%s.tif" % (output_dir, file_name, fov), img_DNAFISH_filtered)

    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH_filtered, blending='additive', colormap='green', contrast_limits=[0, 65535])
    plt.imsave('%s%s_DNAFISH-filtered1_%s.tiff' % (output_dir, file_name, fov), dis.blending(viewer))
    viewer.close()
