import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

mCherry_pos_color = 'gray'
mCherry_neg_color = 'red'
mCherry_pos_threshold = 30000
mCherry_neg_threshold = 20000

sample = 'H2B+BRD1'
file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir2, sample)), na_values=['.'], sep='\t')

n_nuclear_convex_dilation = 2

for fov in range(img_mCherry_stack.shape[0]):
    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    img_mCherry_pos = np.zeros_like(img_seg)
    img_mCherry_neg = np.zeros_like(img_seg)

    mCherry_props = regionprops(img_seg, img_mCherry)
    for i in range(len(mCherry_props)):
        if mCherry_props[i].intensity_mean > mCherry_pos_threshold:
            img_mCherry_pos[img_seg == mCherry_props[i].label] = mCherry_props[i].label
        elif mCherry_props[i].intensity_mean < mCherry_neg_threshold:
            img_mCherry_neg[img_seg == mCherry_props[i].label] = mCherry_props[i].label

    img_mCherry_pos_boundary = dilation(img_mCherry_pos, disk(n_nuclear_convex_dilation+1))
    img_mCherry_pos_boundary[dilation(img_mCherry_pos, disk(n_nuclear_convex_dilation-1)) >= 1] = 0
    img_mCherry_neg_boundary = dilation(img_mCherry_neg, disk(n_nuclear_convex_dilation + 1))
    img_mCherry_neg_boundary[dilation(img_mCherry_neg, disk(n_nuclear_convex_dilation-1)) >= 1] = 0

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bgc.max()])
    viewer.add_image(img_DNAFISH_bgc, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH_bgc.max()])
    viewer.add_image(img_mCherry_pos_boundary, blending='additive', colormap=mCherry_pos_color, contrast_limits=[0, 1])
    viewer.add_image(img_mCherry_neg_boundary, blending='additive', colormap=mCherry_neg_color, contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_%s_img_label.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()