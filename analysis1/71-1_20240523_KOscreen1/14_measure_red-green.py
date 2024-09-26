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
import math
import nd2
import shared.segmentation as seg
import shared.objects as obj
import os
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
total_fov = 16
dshape_factor = 0.145
pixel_size = 300/2720  # uM
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))
img_before_GFP = imutils.rotate(skio.imread("%s/%s/03_%s_GFP_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)
img_before_mCherry = imutils.rotate(skio.imread("%s/%s/03_%s_mCherry_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)

# pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
align = pd.read_csv('%s/%s/08_%s_alignment.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

data = pd.DataFrame()

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    # seg_z = pd_seg_z['seg_z'][fov]
    # img_hoechst_seg_z = img_hoechst[seg_z]
    img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, sample, fov), plugin="tifffile")
    # img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, sample, fov), plugin="tifffile")

    dsize = int(img_hoechst_merge.shape[1] * dshape_factor)
    x = int(align['topleft_x'][fov])
    y = int(align['topleft_y'][fov])
    img_before_GFP_fov = img_before_GFP[y:(y+dsize), x:(x+dsize)]
    img_before_mCherry_fov = img_before_mCherry[y:(y+dsize), x:(x+dsize)]
    img_before_GFP_fov_resize = cv2.resize(img_before_GFP_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)
    img_before_mCherry_fov_resize = cv2.resize(img_before_mCherry_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)

    label_props = regionprops(label(img_seg_total), img_seg_total)
    GFP_props = regionprops(label(img_seg_total), img_before_GFP_fov_resize)
    mCherry_props = regionprops(label(img_seg_total), img_before_mCherry_fov_resize)

    data_temp = pd.DataFrame()
    data_temp['sample'] = [sample] * len(label_props)
    data_temp['fov'] = [fov] * len(label_props)
    data_temp['label_mean_int'] = [label_props[i].intensity_mean for i in range(len(label_props))]
    data_temp['GFP'] = [GFP_props[i].intensity_mean for i in range(len(label_props))]
    data_temp['mCherry'] = [mCherry_props[i].intensity_mean for i in range(len(label_props))]
    data = pd.concat([data, data_temp], axis=0)

    if not os.path.exists("%s/%s/14_red_green_color/" % (output_dir, sample)):
        os.makedirs("%s/%s/14_red_green_color/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
    viewer.add_image(img_before_GFP_fov_resize, blending='additive', colormap='green', contrast_limits=[0, 40000])
    viewer.add_image(img_before_mCherry_fov_resize, blending='additive', colormap='red', contrast_limits=[0, 40000])
    plt.imsave("%s%s/14_red_green_color/%s_%s_red_green_color.tiff" % (output_dir, sample, sample, fov),
               dis.blending(viewer))
    viewer.close()

data.to_csv('%s/%s/14_%s_red_green.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")
