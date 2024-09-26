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
from skimage.measure import label, regionprops
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G11'
samples = ['G11']
total_fovs = [16]

dshape_factor = 0.145
pixel_size = 300/2720  # uM

img_before_GFP = imutils.rotate(skio.imread("%s/%s/03_%s_GFP_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)
img_before_mCherry = imutils.rotate(skio.imread("%s/%s/03_%s_mCherry_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)

align = pd.read_csv('%s/%s/08_%s_alignment.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)
    align_s = align[align['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)
        img_DNAFISH_merge = img_DNAFISH.max(axis=0)
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_hoechst_seg_z = img_hoechst[seg_z]
        img_DNAFISH_seg_z = img_DNAFISH[seg_z]
        img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")

        dsize = int(img_hoechst_seg_z.shape[1] * dshape_factor)
        x = int(align_s['topleft_x'][fov])
        y = int(align_s['topleft_y'][fov])
        img_before_GFP_fov = img_before_GFP[y:(y+dsize), x:(x+dsize)]
        img_before_mCherry_fov = img_before_mCherry[y:(y+dsize), x:(x+dsize)]
        img_before_GFP_fov_resize = cv2.resize(img_before_GFP_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)
        img_before_mCherry_fov_resize = cv2.resize(img_before_mCherry_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)

        label_props = regionprops(label(img_seg_total), img_seg_total)
        GFP_props = regionprops(label(img_seg_total), img_before_GFP_fov_resize)
        mCherry_props = regionprops(label(img_seg_total), img_before_mCherry_fov_resize)

        if not os.path.exists("%s/%s/30_check_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/30_check_color/" % (output_dir, sample))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_seg_z, blending='additive', colormap='blue', contrast_limits=[0, 20000])
        viewer.add_image(img_DNAFISH_seg_z, blending='additive', colormap='green', contrast_limits=[0, 10000])
        plt.imsave("%s%s/30_check_color/30_%s_%s_seg_z.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_before_GFP_fov_resize, blending='additive', colormap='green', contrast_limits=[0, 30000])
        viewer.add_image(img_before_mCherry_fov_resize, blending='additive', colormap='red', contrast_limits=[0, 30000])
        plt.imsave("%s%s/30_check_color/30_%s_%s_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 20000])
        viewer.add_image(img_DNAFISH_merge, blending='additive', colormap='green', contrast_limits=[0, 10000])
        plt.imsave("%s%s/30_check_color/30_%s_%s_merge.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

        """viewer = napari.Viewer()
        viewer.add_image(img_hoechst_seg_z, blending='additive', colormap='blue', contrast_limits=[0, 20000])
        viewer.add_image(img_DNAFISH_seg_z, blending='additive', colormap='green', contrast_limits=[0, 10000])
        viewer.add_image(img_before_GFP_fov_resize, blending='additive', colormap='green', contrast_limits=[0, 30000])
        viewer.add_image(img_before_mCherry_fov_resize, blending='additive', colormap='red', contrast_limits=[0, 30000])
        # plt.imsave("%s%s/30_check_color/30_%s_%s_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        napari.run()"""

print("DONE!")