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
from skimage.filters import threshold_otsu
import shared.segmentation as seg
import shared.objects as obj
import os
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B9'
samples = ['B9']
total_fovs = [16]
dshape_factor = 0.145
pixel_size = 300/2720  # uM
local_size = 200

data = pd.DataFrame(columns=['sample', 'fov', 'DNAFISH_thresh'])
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_hoechst_seg_z = img_hoechst[seg_z]
        img_DNAFISH_seg_z = img_DNAFISH[seg_z]

        img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_ecDNA = np.zeros_like(img_DNAFISH_seg_z)

        thresh = threshold_otsu(img_DNAFISH_seg_z)
        print(thresh)
        img_seg_ecDNA[img_DNAFISH_seg_z > thresh] = 1
        data.loc[len(data.index)] = [s, fov, thresh]

        if not os.path.exists("%s/%s/17_seg_DNAFISH_tif/" % (output_dir, sample)):
            os.makedirs("%s/%s/17_seg_DNAFISH_tif/" % (output_dir, sample))
        tif.imwrite("%s/%s/17_seg_DNAFISH_tif/%s_%s_seg_DNAFISH.tif" % (output_dir, sample, s, fov), img_seg_ecDNA)

        viewer = napari.Viewer()
        viewer.add_image(img_DNAFISH_seg_z, blending='additive', colormap='green', contrast_limits=[0, 20000])
        viewer.add_image(img_seg_ecDNA, blending='additive', contrast_limits=[0, 3])
        if not os.path.exists("%s/%s/17_seg_DNAFISH_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/17_seg_DNAFISH_color/" % (output_dir, sample))
        plt.imsave("%s%s/17_seg_DNAFISH_color/%s_%s_seg_DNAFISH_color.tiff" % (output_dir, sample, s, fov),
                   dis.blending(viewer))
        viewer.close()
        # napari.run()

data.to_csv('%s/%s/17_%s_DNAFISH_thresh.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")
