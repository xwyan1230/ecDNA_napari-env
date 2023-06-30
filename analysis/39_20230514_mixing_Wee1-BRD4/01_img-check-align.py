import skimage.io as skio
import napari
import cv2
import imutils
import shared.display as dis
import shared.image as ima
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '2_Ctrl-GFP-6hr_Ctrl-mCh-6hr'
dshape_factor = 0.47  # 766nm for Keyence 10x, 360nm for 512x512 at 63x
topleft_target = [3333.627431155869, 1430.784821051334]
a_target = 1.3

img_hoechst_Keyence = skio.imread("%s%s/Hoechst.tif" % (data_dir, sample), plugin="tifffile")[:, :, 2]
img_hoechst_DNAFISH = skio.imread("%s%s/%s_%s_hoechst_512_Merging_001_RAW_ch00.tif" % (data_dir, sample, exp, sample), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
# moving down
img_hoechst_DNAFISH2 = np.concatenate([np.zeros(shape=[int(topleft_target[1]), img_hoechst_DNAFISH1.shape[1]]), img_hoechst_DNAFISH1], axis=0)
# moving right
img_hoechst_DNAFISH3 = np.concatenate([np.zeros(shape=[img_hoechst_DNAFISH2.shape[0], int(topleft_target[0])]), img_hoechst_DNAFISH2], axis=1)
img_hoechst_DNAFISH4 = imutils.rotate(img_hoechst_DNAFISH3, angle=a_target)

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH4, blending='additive', colormap='green', contrast_limits=[0, 45535])
# plt.imsave("%s%s/align/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()