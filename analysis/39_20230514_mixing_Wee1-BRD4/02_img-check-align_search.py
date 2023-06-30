import skimage.io as skio
import napari
import cv2
import imutils
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
dshape_factor = 0.47  # 766nm for Keyence 10x, 360nm for 512x512 at 63x
interval = 50

img_hoechst_Keyence = skio.imread("%s%s/Hoechst.tif" % (data_dir, sample), plugin="tifffile")[:, :, 2]
img_hoechst_DNAFISH = skio.imread("%s%s/%s_%s_hoechst_512_Merging_001_RAW_ch00.tif" % (data_dir, sample, exp, sample), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
shapes_layer = viewer.layers['Shapes']
poly_data = shapes.data[0]
topleft_ori = [poly_data[0][1], poly_data[0][0]]
xrange = poly_data[1][1] - poly_data[0][1]
yrange = poly_data[2][0] - poly_data[0][0]
print(xrange)
print(yrange)

# topleft = [25295, 950]
min_ratio = 1
topleft_target = topleft_ori
i = 0
for x in range(int(xrange/interval)):
    for y in range(int(yrange/interval)):
        i = i+1
        topleft = [topleft_ori[0] + x*interval, topleft_ori[1] + y*interval]
        img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
        img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft[1]):int(topleft[1]+img_hoechst_DNAFISH1.shape[0]), int(topleft[0]):int(topleft[0]+img_hoechst_DNAFISH1.shape[1])]

        minus = img_hoechst_DNAFISH1.astype(float) - img_hoechst_Keyence_cut.astype(float)
        minus[minus < 0] = 0
        min_ratio_temp = minus.sum()/img_hoechst_DNAFISH1.sum()
        print('%s/%s: %s' % (i, int(xrange/interval) * int(yrange/interval), min_ratio_temp))
        if min_ratio_temp < min_ratio:
            min_ratio = min_ratio_temp
            topleft_target = [topleft_ori[0] + x*interval, topleft_ori[1] + y*interval]
print(min_ratio)
print(topleft_target)

# topleft_target = [3313.627431155869, 1430.784821051334]  # '2_Ctrl-GFP-6hr_Ctrl-mCh-6hr'

img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft_target[1]):int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH1, blending='additive', colormap='green', contrast_limits=[0, 45535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()

# moving down
img_hoechst_DNAFISH2 = np.concatenate([np.zeros(shape=[int(topleft_target[1]), img_hoechst_DNAFISH1.shape[1]]), img_hoechst_DNAFISH1], axis=0)
# moving right
img_hoechst_DNAFISH3 = np.concatenate([np.zeros(shape=[img_hoechst_DNAFISH2.shape[0], int(topleft_target[0])]), img_hoechst_DNAFISH2], axis=1)
# img_hoechst_DNAFISH4 = imutils.rotate(img_hoechst_DNAFISH3, angle=1.1)

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH3, blending='additive', colormap='green', contrast_limits=[0, 45535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()