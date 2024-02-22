import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.objects as obj
import shared.math as mat
import cv2
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

row = 'C'
sample = 'C6'
total_fov = 16
dshape_factor = 0.0765
n_nuclear_convex_dilation = 4
local_size = 200
rmax = 100
start_fov = 0

df_align = pd.read_csv('%s/alignment/%s/%s_alignment_autocheck.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (data_dir, sample, sample), plugin="tifffile")

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear_bgc = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    s = img_nuclear_bgc.shape[0]*dshape_factor
    topleft = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]
    img_before_GFP_cut = img_before_GFP.copy()[int(topleft[1]):int(topleft[1]+s)+5, int(topleft[0]):int(topleft[0]+s)+5]
    s1 = img_before_GFP_cut.shape[0]
    img_before_GFP_cut_resize = cv2.resize(img_before_GFP_cut, dsize=(int(s1 * 1/dshape_factor), int(s1 * 1/dshape_factor)),
                                           interpolation=cv2.INTER_AREA)[:img_nuclear_bgc.shape[0], :img_nuclear_bgc.shape[1]]
    img_before_mCherry_cut = img_before_mCherry.copy()[int(topleft[1]):int(topleft[1]+s)+5, int(topleft[0]):int(topleft[0]+s)+5]
    img_before_mCherry_cut_resize = cv2.resize(img_before_mCherry_cut,
                                           dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                           interpolation=cv2.INTER_AREA)[:img_nuclear_bgc.shape[0], :img_nuclear_bgc.shape[1]]

    viewer = napari.Viewer()
    viewer.add_image(img_before_mCherry_cut_resize, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_before_GFP_cut_resize, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH_bgc, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()

print("DONE!")