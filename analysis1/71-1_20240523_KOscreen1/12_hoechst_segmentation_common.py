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
data = pd.read_csv("%s/%s/11_%s_seg_centroids.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')
data_common = pd.DataFrame(columns=['sample', 'fov', 'label', 'seg', 'centroid_0', 'centroid_1', 'mean_int', 'label_new'])

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    data_total = data[(data['fov'] == fov) & (data['seg'] == 'total')].copy().reset_index(drop=True)
    data_z = data[(data['fov'] == fov) & (data['seg'] == 'z')].copy().reset_index(drop=True)
    count = 1
    for i in range(len(data_total)):
        print("%s/%s" % (i+1, len(data_total)))
        for j in range(len(data_z)):
            if (data_z['centroid_0'][j]>(data_total['centroid_0'][i]-5)) & (data_z['centroid_0'][j]<(data_total['centroid_0'][i]+5)) & (data_z['centroid_1'][j]>(data_total['centroid_1'][i]-5)) & (data_z['centroid_1'][j]<(data_total['centroid_1'][i]+5)):
                data_common.loc[len(data_common.index)] = [sample, fov, data_total['label'][i], 'total', data_total['centroid_0'][i], data_total['centroid_1'][i], data_total['mean_int'][i], count]
                data_common.loc[len(data_common.index)] = [sample, fov, data_z['label'][j], 'z', data_z['centroid_0'][j], data_z['centroid_1'][j], data_z['mean_int'][j], count]
                count = count + 1
                break

data_common.to_csv('%s/%s/12_%s_seg_common_centroids.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")