import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM'
dshape_factor = 0.392  # 142nm for 2048x2048 40x, 360nm for 512x512 at 63x, by calculation this factor is roughly 0.394

# left
img_hoechst_IF = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
ori_shape = img_hoechst_IF.shape
img_hoechst_IF_resize = cv2.resize(img_hoechst_IF, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
img_hoechst_IF_resize1 = np.concatenate([np.zeros(shape=[820, img_hoechst_IF_resize.shape[1]]), img_hoechst_IF_resize], axis=0)
img_hoechst_DNAFISH = skio.imread("%s%s/DM_hoechst_DNAFISH_merge.tif" % (data_dir, sample), plugin="tifffile")
img_hoechst_DNAFISH1 = np.concatenate([np.zeros(shape=[img_hoechst_DNAFISH.shape[0], 700]), img_hoechst_DNAFISH], axis=1)


viewer = napari.Viewer()
viewer.add_image(img_hoechst_IF_resize1, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH1, blending='additive', colormap='green', contrast_limits=[0, 45535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()