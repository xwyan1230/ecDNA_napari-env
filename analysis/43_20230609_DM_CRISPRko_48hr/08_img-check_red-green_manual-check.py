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
sample = 'C2'
total_fov = 16
dshape_factor = 0.0765
n_col = 4
border = 50
interval = 1

df_align = pd.read_csv('%s/alignment/%s_alignment.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y'])

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_hoechst_DNAFISH_cut = skio.imread("%s/FISH/%s_hoechst_DNAFISH_withborder.tif" % (data_dir, sample), plugin='tifffile')

for fov in range(total_fov):
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_hoechst = skio.imread("%s/FISH/%s_ch01.tif" % (data_dir, filename), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    topleft = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]

    count = 0
    while count == 0:
        print(topleft)
        img_hoechst_local = np.zeros_like(img_after_hoechst)
        img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize,[int(topleft[1]), int(topleft[0])])
        viewer = napari.Viewer()
        viewer.add_image(img_before_GFP, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_before_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
        # viewer.add_image(img_after_hoechst, blending='additive', contrast_limits=[0, 65535])
        viewer.add_image(img_hoechst_local, blending='additive', colormap='blue', contrast_limits=[0, 35535])
        shapes_up = viewer.add_shapes(name='up', ndim=2)
        shapes_down = viewer.add_shapes(name='down', ndim=2)
        shapes_left = viewer.add_shapes(name='left', ndim=2)
        shapes_right = viewer.add_shapes(name='right', ndim=2)
        napari.run()

        if len(shapes_right.data) + len(shapes_left.data) + len(shapes_up.data) + len(shapes_down.data) == 0:
            count = 1
        else:
            topleft_modify = [topleft[0]+interval*(len(shapes_right.data)-len(shapes_left.data)), topleft[1]+interval*(len(shapes_down.data)-len(shapes_up.data))]
            topleft = topleft_modify
    data.loc[len(data.index)] = [sample, fov, topleft[0], topleft[1]]
    data.to_csv('%s/%s_alignment_check.txt' % (output_dir, sample), index=False, sep='\t')
    print(topleft)

print("DONE!")
