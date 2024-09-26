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

data = pd.DataFrame(columns=['sample', 'fov', 'label', 'seg', 'centroid_0', 'centroid_1', 'mean_int'])

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_seg_total = skio.imread("%s/%s/09_seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov), plugin="tifffile")
    img_seg_z = skio.imread("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, sample, fov), plugin="tifffile")

    total_props = regionprops(label(img_seg_total), img_seg_total)
    for i in range(len(total_props)):
        data.loc[len(data.index)] = [sample, fov, total_props[i].label, 'total', total_props[i].centroid[0], total_props[i].centroid[1], total_props[i].intensity_mean]
    z_props = regionprops(label(img_seg_z), img_seg_z)
    for i in range(len(z_props)):
        data.loc[len(data.index)] = [sample, fov, z_props[i].label, 'z', z_props[i].centroid[0], z_props[i].centroid[1], z_props[i].intensity_mean]

data.to_csv('%s/%s/11_%s_seg_centroids.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")