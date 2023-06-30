import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import math
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
from skimage.measure import label, regionprops_table, regionprops
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os
import skimage.io as skio

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'

img_hoechst_Keyence = skio.imread("%s%s/Hoechst_cut.tif" % (data_dir, sample), plugin="tifffile")
img_nuclear_seg = skio.imread("%s%s/seg_keyence/Hoechst_keyence_seg_convex.tif" % (data_dir1, sample), plugin="tifffile")
img_mCherry_Keyence = skio.imread("%s%s/mCherry_cut.tif" % (data_dir, sample), plugin="tifffile")
img_GFP_Keyence = skio.imread("%s%s/GFP_cut.tif" % (data_dir, sample), plugin="tifffile")

nuclear_props_hoechst = regionprops(img_nuclear_seg, img_hoechst_Keyence)
nuclear_props_mCherry = regionprops(img_nuclear_seg, img_mCherry_Keyence)
nuclear_props_GFP = regionprops(img_nuclear_seg, img_GFP_Keyence)

data = pd.DataFrame(columns=['nuclear', 'label', 'area_nuclear', 'mean_int_nuclear', 'mean_int_GFP', 'mean_int_mCherry'])
data['nuclear'] = range(len(nuclear_props_mCherry))
data['label'] = [nuclear_props_mCherry[x].label for x in range(len(nuclear_props_mCherry))]
data['area_nuclear'] = [nuclear_props_mCherry[x].area for x in range(len(nuclear_props_mCherry))]
data['mean_int_nuclear'] = [nuclear_props_hoechst[x].intensity_mean for x in range(len(nuclear_props_hoechst))]
data['mean_int_GFP'] = [nuclear_props_GFP[x].intensity_mean for x in range(len(nuclear_props_GFP))]
data['mean_int_mCherry'] = [nuclear_props_mCherry[x].intensity_mean for x in range(len(nuclear_props_mCherry))]

data.to_csv('%s%s/analysis_keyence.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")