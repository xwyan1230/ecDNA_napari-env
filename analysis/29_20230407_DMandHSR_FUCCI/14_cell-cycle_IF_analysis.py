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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'


def image_stdev(region, intensities):
    # note the ddof arg to get the sample var if you so desire!
    return np.std(intensities[region], ddof=1)


img_hoechst = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_nuclear_seg = skio.imread("%s%s/DM_nuclear_red_green_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")
img_red = skio.imread("%s%s/DM_red_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_green = skio.imread("%s%s/DM_green_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_MYC = skio.imread("%s%s/DM_MYC_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")

nuclear_props_hoechst = regionprops(img_nuclear_seg, img_hoechst, extra_properties=[image_stdev])
nuclear_props_red = regionprops(img_nuclear_seg, img_red)
nuclear_props_green = regionprops(img_nuclear_seg, img_green)
nuclear_props_MYC = regionprops(img_nuclear_seg, img_MYC)

data = pd.DataFrame(columns=['nuclear', 'label', 'area_nuclear', 'mean_int_nuclear', 'mean_int_green', 'mean_int_red', 'mean_int_MYC', 'intensity_stdev'])
data['nuclear'] = range(len(nuclear_props_red))
data['label'] = [nuclear_props_red[x].label for x in range(len(nuclear_props_red))]
data['area_nuclear'] = [nuclear_props_red[x].area for x in range(len(nuclear_props_red))]
data['mean_int_nuclear'] = [nuclear_props_hoechst[x].intensity_mean for x in range(len(nuclear_props_hoechst))]
data['mean_int_green'] = [nuclear_props_green[x].intensity_mean for x in range(len(nuclear_props_green))]
data['mean_int_red'] = [nuclear_props_red[x].intensity_mean for x in range(len(nuclear_props_red))]
data['mean_int_MYC'] = [nuclear_props_MYC[x].intensity_mean for x in range(len(nuclear_props_MYC))]
data['intensity_stdev'] = [nuclear_props_hoechst[x].image_stdev for x in range(len(nuclear_props_hoechst))]

data.to_csv('%s%s_cellcycle.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")