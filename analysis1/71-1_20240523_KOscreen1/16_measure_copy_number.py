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

data = pd.DataFrame()
df_bg = pd.read_csv('%s/bg.txt' % output_dir, na_values=['.'], sep='\t')

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    img_DNAFISH = img_stack[fov, :, 1, :, :]
    img_DNAFISH_sum = img_DNAFISH.astype(int).sum(axis=0)

    img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, sample, fov), plugin="tifffile")
    label_props = regionprops(label(img_seg_total), img_seg_total)
    DNAFISH_props = regionprops(label(img_seg_total), img_DNAFISH_sum)
    bg_merge = df_bg[df_bg['sample'] == sample]['bg_merge'].tolist()[0]

    data_temp = pd.DataFrame()
    data_temp['sample'] = [sample] * len(label_props)
    data_temp['fov'] = [fov] * len(label_props)
    data_temp['label_mean_int'] = [label_props[i].intensity_mean for i in range(len(label_props))]
    data_temp['nuclear_area_merge'] = [label_props[i].area for i in range(len(label_props))]
    data_temp['DNAFISH_mean_int_merge'] = [DNAFISH_props[i].intensity_mean for i in range(len(DNAFISH_props))]
    data_temp['bg_merge'] = [bg_merge] * len(label_props)

    data = pd.concat([data, data_temp], axis=0)

data['DNAFISH_total_int_merge'] = [(data['DNAFISH_mean_int_merge'].tolist()[i]-data['bg_merge'].tolist()[i]) *
                                   data['nuclear_area_merge'].tolist()[i] if data['DNAFISH_mean_int_merge'].tolist()[i]-data['bg_merge'].tolist()[i] >0
                                   else 0 for i in range(len(data))]

data.to_csv('%s/%s/16_%s_copy_number.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")
