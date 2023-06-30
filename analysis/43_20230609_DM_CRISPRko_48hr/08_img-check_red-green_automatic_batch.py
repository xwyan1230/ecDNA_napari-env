import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%salignment/" % master_folder

row = 'C'
sample = 'C3'
batch = 2
total_fov = 12
dshape_factor = 0.0765

df_align = pd.read_csv('%s/alignment/%s/%s_alignment_%s.txt' % (master_folder, sample, sample, batch), na_values=['.'], sep='\t')
data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_hoechst_DNAFISH_cut = skio.imread("%s/FISH/%s/%s_hoechst_DNAFISH_withborder.tif" % (data_dir, sample, sample), plugin='tifffile')
img_before = img_before_GFP.copy()
img_before[img_before_mCherry>img_before_GFP] = img_before_mCherry[img_before_mCherry>img_before_GFP]
img_hoechst_local = np.zeros_like(img_after_hoechst)

for fov in range(total_fov):
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s%s' % (row, sample, fov)
    img_hoechst = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    topleft = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]
    s = img_hoechst_resize.shape[0]
    print(topleft)
    topleft_target, min_ratio = ima.img_search_local(img_before, img_hoechst_resize, topleft, 1)
    data.loc[len(data.index)] = [sample, fov, topleft_target[0], topleft_target[1], min_ratio]
    print(min_ratio)
    print(topleft_target)

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize,
                                           [int(topleft_target[1]), int(topleft_target[0])])

    data.to_csv('%s/%s/%s_alignment_autocheck_%s.txt' % (output_dir, sample, sample, batch), index=False, sep='\t')
    viewer = napari.Viewer()
    viewer.add_image(img_before, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()

viewer = napari.Viewer()
viewer.add_image(img_before, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave("%s/%s/%s_alignment_local_before_%s.tiff" % (output_dir, sample, sample, batch), dis.blending(viewer))
viewer.close()

print("DONE!")
