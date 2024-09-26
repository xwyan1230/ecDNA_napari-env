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
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))

data = pd.read_csv("%s/%s/12_%s_seg_common_centroids.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for fov in range(total_fov):
    print("%s/%s" % (fov+1, total_fov))
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    seg_z = pd_seg_z['seg_z'][fov]
    img_hoechst_seg_z = img_hoechst[seg_z]
    img_seg_total = skio.imread("%s/%s/09_seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov), plugin="tifffile")
    img_seg_z = skio.imread("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, sample, fov), plugin="tifffile")
    img_seg_total_new = np.zeros_like(img_seg_total)
    img_seg_z_new = np.zeros_like(img_seg_z)
    data_total = data[(data['fov'] == fov) & (data['seg'] == 'total')].copy().reset_index(drop=True)
    data_z = data[(data['fov'] == fov) & (data['seg'] == 'z')].copy().reset_index(drop=True)

    total_props = regionprops(label(img_seg_total), img_seg_total)
    for i in range(len(total_props)):
        if total_props[i].intensity_mean in data_total['mean_int'].tolist():
            label_index = data_total[data_total['mean_int'] == total_props[i].intensity_mean].index
            print('total: %s-%s' % (total_props[i].intensity_mean, data_total['label_new'][label_index]))
            img_seg_total_new[img_seg_total == total_props[i].intensity_mean] = data_total['label_new'][label_index]

    z_props = regionprops(label(img_seg_z), img_seg_z)
    for i in range(len(z_props)):
        if z_props[i].intensity_mean in data_z['mean_int'].tolist():
            label_index = data_z[data_z['mean_int'] == z_props[i].intensity_mean].index
            print('z: %s-%s' % (z_props[i].intensity_mean, data_z['label_new'][label_index]))
            img_seg_z_new[img_seg_z == z_props[i].intensity_mean] = data_z['label_new'][label_index]

    """test = img_seg_total_new-img_seg_z_new
    viewer = napari.Viewer()
    viewer.add_image(img_seg_total_new, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_seg_z_new, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer.add_image(test, blending='additive', contrast_limits=[0, 1])
    napari.run()"""

    if not os.path.exists("%s/%s/13_seg_new_tif/" % (output_dir, sample)):
        os.makedirs("%s/%s/13_seg_new_tif/" % (output_dir, sample))
    tif.imwrite("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, sample, fov), img_seg_total_new)
    tif.imwrite("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, sample, fov), img_seg_z_new)

    if not os.path.exists("%s/%s/13_seg_new_color/" % (output_dir, sample)):
        os.makedirs("%s/%s/13_seg_new_color/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
    viewer.add_image(img_seg_total_new, blending='additive', contrast_limits=[0, 1])
    plt.imsave("%s%s/13_seg_new_color/%s_%s_seg_new_color.tiff" % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_seg_z, blending='additive', colormap='blue', contrast_limits=[0, 12000])
    viewer.add_image(img_seg_z_new, blending='additive', contrast_limits=[0, 1])
    plt.imsave("%s%s/13_seg_new_color/%s_%s_seg_z_new_color.tiff" % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")