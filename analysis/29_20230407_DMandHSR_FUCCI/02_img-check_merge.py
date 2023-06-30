import skimage.io as skio
import napari
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

sample = 'DM_DNAFISH_hoechst'

# left
img_hoechst_NW = skio.imread("%s%s/DM_hoechst_NW_merge_RAW_ch00.tif" % (data_dir, sample), plugin="tifffile")
# print(img_hoechst_NW.shape)
img_hoechst_NWSW = skio.imread("%s%s/DM_hoechst_NW-SW_merge_RAW_ch00.tif" % (data_dir, sample), plugin="tifffile")
# print(img_hoechst_NWSW.shape)
img_hoechst1 = np.concatenate([img_hoechst_NW[:(4626-536), :(6043-215)], img_hoechst_NWSW[:1414, 215:6043]], axis=0)
img_hoechst_check1 = img_hoechst_NW
# print(img_hoechst1.shape)
img_hoechst_SW = skio.imread("%s%s/DM_hoechst_SW_merge_RAW_ch00.tif" % (data_dir, sample), plugin="tifffile")
# print(img_hoechst_SW.shape)
img_hoechst2 = np.concatenate([img_hoechst1[:5504, :(5828-25)], img_hoechst_SW[615:, 25:5828]], axis=0)
img_hoechst_check2 = np.concatenate([np.zeros(shape=[5504-615, 5828-25]), img_hoechst_SW[:, 25:5828]], axis=0)
# print(img_hoechst2.shape())

# right
img_hoechst_NE = skio.imread("%s%s/DM_hoechst_NE_merge_RAW_ch00.tif" % (data_dir, sample), plugin="tifffile")
# print(img_hoechst_NE.shape)
img_hoechst_SE = skio.imread("%s%s/DM_hoechst_SE_merge_RAW_ch00.tif" % (data_dir, sample), plugin="tifffile")
# print(img_hoechst_SE.shape)
img_hoechst3 = np.concatenate([img_hoechst_NE[:(4635-213), :(6050-28)], img_hoechst_SE[:4649, 28:6050]], axis=0)
img_hoechst_check3 = img_hoechst_NE

# combine
# print(img_hoechst2.shape)
# print(img_hoechst3.shape)
img_hoechst4 = np.concatenate([img_hoechst2[72:(9071+72), :(5803-123)], img_hoechst3[:9071, 20:6022]], axis=1)
img_hoechst_check4 = img_hoechst2[72:, :]
print(img_hoechst4.shape)
tif.imwrite("%s%s/DM_hoechst_merge.tif" % (output_dir, sample), img_hoechst4)

viewer = napari.Viewer()
viewer.add_image(img_hoechst4, blending='additive', colormap='blue', contrast_limits=[0, 35535])
# viewer.add_image(img_hoechst_check4, blending='additive', colormap='green', contrast_limits=[0, 45535])
plt.imsave("%s%s/DM_hoechst_merge.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()


