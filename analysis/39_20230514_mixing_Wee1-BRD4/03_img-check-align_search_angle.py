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
sample = '6_Wee1-GFP-6hr_Ctrl-mCh-6hr'
dshape_factor = 0.47  # 766nm for Keyence 10x, 360nm for 512x512 at 63x
interval = 10

# topleft_target = [3313.627431155869, 1430.784821051334]  2_Ctrl-GFP-6hr_Ctrl-mCh-6hr [3333.627431155869, 1430.784821051334] 1.3
# topleft_target = [3753.8435086956088, 5124.393828360026]  5_Ctrl-GFP-6hr_Wee1-mCh-6hr [3733.8435086956088, 5094.393828360026] 1.1
# topleft_target = [5187.786656194921, 4282.431129459471]  6_Wee1-GFP-6hr_Ctrl-mCh-6hr [5177.786656194921, 4282.431129459471] 1.1
# topleft_target = [2615.5115503413313, 4631.470673742484]  9_Ctrl-GFP-6hr_BRD4-mCh-6hr [2635.5115503413313, 4641.470673742484] 0.5
# topleft_target = [3875.6283976329314, 2309.720277882853]  10_BRD4-GFP-6hr_Ctrl-mCh-6hr [3865.6283976329314, 2319.720277882853] 0.3

img_hoechst_Keyence = skio.imread("%s%s/Hoechst.tif" % (data_dir, sample), plugin="tifffile")[:, :, 2]
img_hoechst_DNAFISH = skio.imread("%s%s/%s_%s_hoechst_512_Merging_001_RAW_ch00.tif" % (data_dir, sample, exp, sample), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)

xrange = 100
yrange = 100
arange = 2
a_interval = 0.1
min_ratio = 1
topleft_target = [5187.786656194921, 4282.431129459471]
topleft_ori = [topleft_target[0]-xrange/2, topleft_target[1]-yrange/2]
a_target = 0
i = 0
for x in range(int(xrange/interval)):
    for y in range(int(yrange/interval)):
        for a in range(int(arange/a_interval)):
            i = i+1
            topleft = [topleft_ori[0] + x*interval, topleft_ori[1] + y*interval]
            img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
            img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft[1]):int(topleft[1]+img_hoechst_DNAFISH1.shape[0]), int(topleft[0]):int(topleft[0]+img_hoechst_DNAFISH1.shape[1])]
            img_hoechst_DNAFISH_rotate = imutils.rotate(img_hoechst_DNAFISH1, angle=a*a_interval)

            minus = img_hoechst_DNAFISH_rotate.astype(float) - img_hoechst_Keyence_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum()/img_hoechst_DNAFISH_rotate.sum()
            print('%s/%s: %s' % (i, int(xrange/interval) * int(yrange/interval) * int(arange/a_interval), min_ratio_temp))
            if min_ratio_temp < min_ratio:
                min_ratio = min_ratio_temp
                topleft_target = [topleft_ori[0] + x*interval, topleft_ori[1] + y*interval]
                a_target = a*a_interval
print(sample)
print(min_ratio)
print(topleft_target)
print(a_target)

img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft_target[1]):int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0]), int(topleft_target[0]):int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])]
img_hoechst_DNAFISH_rotate = imutils.rotate(img_hoechst_DNAFISH1, angle=a_target)

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH_rotate, blending='additive', colormap='green', contrast_limits=[0, 45535])
if not os.path.exists("%s%s/align/" % (output_dir, sample)):
    os.makedirs("%s%s/align/" % (output_dir, sample))
plt.imsave("%s%s/align/DM_alignment_local.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()

# moving down
img_hoechst_DNAFISH2 = imutils.rotate(img_hoechst_DNAFISH1, angle=a_target)
img_hoechst_DNAFISH3 = ima.image_paste(img_hoechst_Keyence, img_hoechst_DNAFISH2, [int(topleft_target[1]), int(topleft_target[0])])

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_DNAFISH3, blending='additive', colormap='green', contrast_limits=[0, 45535])
plt.imsave("%s%s/align/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()