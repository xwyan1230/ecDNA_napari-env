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
output_dir = "%sdata/" % master_folder

row = 'C'
sample = 'C2'
total_fov = 16
dshape_factor = 0.0765
n_col = 4
border = 50

df_align = pd.read_csv('%s/alignment/%s_alignment.txt' % (master_folder, sample), na_values=['.'], sep='\t')

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
    img_DNAFISH = skio.imread("%s/FISH/%s_ch00.tif" % (data_dir, filename), plugin="tifffile")
    s = img_hoechst.shape[0]*dshape_factor
    topleft = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]
    img_before_GFP_cut = img_before_GFP.copy()[int(topleft[1])-border:int(topleft[1]+s)+border, int(topleft[0])-border:int(topleft[0]+s)+border]
    s1 = img_before_GFP_cut.shape[0]
    img_before_GFP_cut_resize = cv2.resize(img_before_GFP_cut, dsize=(int(s1 * 1/dshape_factor), int(s1 * 1/dshape_factor)),
                                           interpolation=cv2.INTER_AREA)
    img_after_hoechst_cut = img_after_hoechst.copy()[int(topleft[1])-border:int(topleft[1]+s)+border, int(topleft[0])-border:int(topleft[0]+s)+border]
    img_after_hoechst_cut_resize = cv2.resize(img_after_hoechst_cut, dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                           interpolation=cv2.INTER_AREA)
    img_before_mCherry_cut = img_before_mCherry.copy()[int(topleft[1])-border:int(topleft[1]+s)+border, int(topleft[0])-border:int(topleft[0]+s)+border]
    img_before_mCherry_cut_resize = cv2.resize(img_before_mCherry_cut,
                                           dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                           interpolation=cv2.INTER_AREA)
    img_hoechst_border = ima.image_paste(img_after_hoechst_cut_resize, img_hoechst, [int(border* 1 / dshape_factor), int(border* 1 / dshape_factor)])
    img_DNAFISH_border = ima.image_paste(img_after_hoechst_cut_resize, img_DNAFISH,
                                         [int(border * 1 / dshape_factor), int(border * 1 / dshape_factor)])
    img_red = np.zeros_like(img_hoechst_border)
    img_green = np.zeros_like(img_hoechst_border)

    viewer = napari.Viewer()
    viewer.add_image(img_before_GFP_cut_resize, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_before_mCherry_cut_resize, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_after_hoechst_cut_resize, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_border, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH_border, blending='additive', colormap='green', contrast_limits=[0, 65535])
    shapes_red = viewer.add_shapes(name='red', ndim=2)
    shapes_green = viewer.add_shapes(name='green', ndim=2)
    napari.run()

    """img_red = ima.napari_add_or_remove(shapes_red.data, 'add', img_red)
    img_green = ima.napari_add_or_remove(shapes_green.data, 'add', img_green)
    s2 = img_hoechst.shape[0]
    img_red = img_red[int(border* 1 / dshape_factor): int(border* 1 / dshape_factor)+s2, int(border* 1 / dshape_factor): int(border* 1 / dshape_factor)+s2]
    img_green = img_green[int(border* 1 / dshape_factor): int(border* 1 / dshape_factor)+s2, int(border* 1 / dshape_factor): int(border* 1 / dshape_factor)+s2]
    if not os.path.exists("%s/red-green/%s/" % (output_dir, sample)):
        os.makedirs("%s/red-green/%s/" % (output_dir, sample))
    tif.imwrite("%s/red-green/%s/%s_%s_GFP.tif" % (output_dir, sample, fov, sample), img_green)
    tif.imwrite("%s/red-green/%s/%s_%s_mCherry.tif" % (output_dir, sample, fov, sample), img_red)

    viewer = napari.Viewer()
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 1])
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()"""