import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir1 = "%sdata/singleZ/" % master_folder
data_dir2 = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

mCherry_pos_color = 'red'
mCherry_neg_color = 'gray'
na_color = 'yellow'
mCherry_pos_threshold = 20000
mCherry_neg_threshold = 11000
n_nuclear_convex_dilation = 4
start_fov = 0

sample = 'DM-Ctrl_mix_mCh-BRD3'
batch = 2
if batch > 0:
    file_name = '%s_%s_RAW' % (sample, batch)
else:
    file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir1, sample, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch02.tif" % (data_dir1, sample,  file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir1, sample, file_name), plugin="tifffile")

df = pd.read_csv(("%s%s_n%s.txt" % (data_dir2, sample, n_nuclear_convex_dilation)), na_values=['.'], sep='\t')
max_area = 60000
min_circ = 0.8

for f in range(img_mCherry_stack.shape[0]):
    fov = f + start_fov
    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry_bgc = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir2, sample, batch, sample, fov), plugin="tifffile")

    img_mCherry_pos = np.zeros_like(img_seg)
    img_mCherry_neg = np.zeros_like(img_seg)
    img_mCherry_na = np.zeros_like(img_seg)

    mCherry_props = regionprops(img_seg, img_mCherry_bgc)
    for i in range(len(mCherry_props)):
        if (mCherry_props[i].area < max_area) & ((4 * math.pi * mCherry_props[i].area) / (mCherry_props[i].perimeter ** 2) > min_circ):
            if mCherry_props[i].intensity_mean > mCherry_pos_threshold:
                img_mCherry_pos[img_seg == mCherry_props[i].label] = mCherry_props[i].label
            elif mCherry_props[i].intensity_mean < mCherry_neg_threshold:
                img_mCherry_neg[img_seg == mCherry_props[i].label] = mCherry_props[i].label
            else:
                img_mCherry_na[img_seg == mCherry_props[i].label] = mCherry_props[i].label

    img_mCherry_pos_boundary = dilation(img_mCherry_pos, disk(n_nuclear_convex_dilation+1))
    img_mCherry_pos_boundary[dilation(img_mCherry_pos, disk(n_nuclear_convex_dilation-2)) >= 1] = 0
    img_mCherry_neg_boundary = dilation(img_mCherry_neg, disk(n_nuclear_convex_dilation + 1))
    img_mCherry_neg_boundary[dilation(img_mCherry_neg, disk(n_nuclear_convex_dilation-2)) >= 1] = 0
    img_mCherry_na_boundary = dilation(img_mCherry_na, disk(n_nuclear_convex_dilation + 1))
    img_mCherry_na_boundary[dilation(img_mCherry_na, disk(n_nuclear_convex_dilation - 2)) >= 1] = 0

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bgc.max()])
    viewer.add_image(img_mCherry_bgc, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH_bgc, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH_bgc.max()])
    viewer.add_image(img_mCherry_pos_boundary, blending='additive', colormap=mCherry_pos_color, contrast_limits=[0, 1])
    viewer.add_image(img_mCherry_neg_boundary, blending='additive', colormap=mCherry_neg_color, contrast_limits=[0, 1])
    viewer.add_image(img_mCherry_na_boundary, blending='additive', colormap=na_color, contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s/%s_%s_img_label_withmCh.tiff' % (output_dir, sample, batch, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bgc.max()])
    viewer.add_image(img_DNAFISH_bgc, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH_bgc.max()])
    viewer.add_image(img_mCherry_pos_boundary, blending='additive', colormap=mCherry_pos_color, contrast_limits=[0, 1])
    viewer.add_image(img_mCherry_neg_boundary, blending='additive', colormap=mCherry_neg_color, contrast_limits=[0, 1])
    viewer.add_image(img_mCherry_na_boundary, blending='additive', colormap=na_color, contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s/%s_%s_img_label.tiff' % (output_dir, sample, batch, sample, fov), dis.blending(viewer))
    viewer.close()