import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis, erosion
import shared.objects as obj
import math
import shared.image as ima
import shared.dataframe as dat
from skimage import segmentation
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
import seaborn as sns
import skimage.io as skio
import shared.math as mat
import shared.display as dis
import pandas as pd
import tifffile as tif
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/E/"
sample_lst = ['HZ', 'FM', 'KW', 'FD', 'ME', 'CK', 'AO', 'IK', 'CR', 'BD', 'HD', 'NP', 'BQ', 'CO', 'LM', 'GP',
              'MM', 'LY', 'JZ', 'NB', 'LT', 'MK', 'DA', 'LX', 'GR', 'NQ', 'JF', 'DT']
save_folder = master_folder

for sample in sample_lst:
    # load images
    img_nuclear = skio.imread("%s%s_b.BMP" % (master_folder, sample))[:, :, 2]
    img_DNAFISH = skio.imread("%s%s_g.BMP" % (master_folder, sample))[:, :, 1]
    img_nuclear_seg = skio.imread("%s%s_seg.tif" % (master_folder, sample), plugin="tifffile")
    img_DNAFISH_seg = skio.imread("%s%s_ecSeg.tif" % (master_folder, sample), plugin="tifffile")

    # bg_correction
    _, sample_img_nuclear = seg.get_bg_img(img_nuclear)
    _, sample_img_DNAFISH = seg.get_bg_img(img_DNAFISH)

    sample_seg = np.zeros_like(img_nuclear)
    sample_seg[sample_img_nuclear == 1] = 1
    sample_seg[sample_img_DNAFISH == 1] = 1
    sample_seg = binary_dilation(sample_seg, disk(50))
    bg = np.ones_like(sample_seg)
    bg[sample_seg == 1] = 0

    bg_int_nuclear = np.sum(bg * img_nuclear) / np.sum(bg)
    bg_int_DNAFISH = np.sum(bg * img_DNAFISH) / np.sum(bg)
    img_nuclear_bgc = img_nuclear
    img_DNAFISH_bgc = img_DNAFISH

    img_nuclear_bgc = img_nuclear.astype(float) - np.ones_like(img_nuclear) * bg_int_nuclear
    img_nuclear_bgc[img_nuclear_bgc < 0] = 0
    img_DNAFISH_bgc = img_DNAFISH.astype(float) - np.ones_like(img_DNAFISH) * bg_int_DNAFISH
    img_DNAFISH_bgc[img_DNAFISH_bgc < 0] = 0

    img_nuclear_center = erosion(img_nuclear_seg)
    img_nuclear_border = img_nuclear_seg.copy()
    img_nuclear_border[img_nuclear_center > 0] = 0

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue')
    viewer.add_image(img_DNAFISH_bgc, blending='additive', colormap='green')
    viewer.add_image(img_nuclear_border, blending='additive', contrast_limits=[0, 1])
    plt.imsave('%s%s_analyze.tiff' % (save_folder, sample), dis.blending(viewer))
    viewer.close()

print("DONE!")