import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk, binary_erosion, dilation, erosion
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240213_analysis_GFPandmCherry_IFandFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
threshold = 0
circ_threshold = 0.5
min_size = 300
max_size = 1500

folder = 'RNAFISH_96well'
sample = 'DM_5000'
sample_label = 'XY01'
n_img = 9

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
GFP_lst = []
mCherry_lst = []
RNAFISH_lst = []


def img_crop(img, i):
    if i == 0:
        out = img[0:1440, 0:1920]
    elif (i == 1) | (i == 2):
        out = img[0:1440, 576:1920]
    elif (i == 3) | (i == 6):
        out = img[432:1440, 0:1920]
    elif (i == 4) | (i == 5):
        out = img[432:1440, 0:1344]
    elif (i == 7) | (i == 8):
        out = img[432:1440, 576:1920]
    return out


for i in range(n_img):
    print(i)
    file_name = 'Image_%s_0000%s' % (sample_label, i+1)
    img_hoechst = img_crop(skio.imread("%s%s/%s/%s_CH1.tif" % (data_dir, folder, sample_label, file_name), plugin="tifffile")[:, :, 2], i)
    img_green = img_crop(skio.imread("%s%s/%s/%s_CH2.tif" % (data_dir, folder, sample_label, file_name), plugin="tifffile")[:, :, 1], i)
    img_red = img_crop(skio.imread("%s%s/%s/%s_CH3.tif" % (data_dir, folder, sample_label, file_name), plugin="tifffile")[:, :, 0], i)
    img_farred = img_crop(skio.imread("%s%s/%s/%s_CH4.tif" % (data_dir, folder, sample_label, file_name), plugin="tifffile")[:, :, 0], i)
    img_nuclear_seg = seg.cell_seg_fluorescent(img_hoechst, otsu_factor=otsu_factor, max_size=max_size, maxima_threshold=5,
                                               min_size=min_size, circ_thresh=circ_threshold, threshold=threshold).astype(int)
    print(np.max(img_nuclear_seg))
    img_nuclear_seg = erosion(img_nuclear_seg)
    nuclear_props = regionprops(img_nuclear_seg)
    RNA_mask = dilation(img_nuclear_seg, disk(5))
    RNA_mask[img_nuclear_seg >= 1] = 0

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    viewer.add_image(img_nuclear_seg, blending='additive', colormap='red', contrast_limits=[0, np.max(img_nuclear_seg)])
    viewer.add_image(RNA_mask, blending='additive', colormap='green', contrast_limits=[0, np.max(RNA_mask)])
    napari.run()

    """for j in range(np.max(img_nuclear_seg)):
        label_img_temp = np.zeros_like(img_nuclear_seg)
        label_img_temp[img_nuclear_seg == j] = 1
        RNA_mask_temp = np.zeros_like(img_nuclear_seg)
        RNA_mask_temp[RNA_mask == j] = 1
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
        viewer.add_image(img_farred, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        viewer.add_image(label_img_temp, blending='additive', colormap='red', contrast_limits=[0, 1])
        viewer.add_image(RNA_mask_temp, blending='additive', colormap='green', contrast_limits=[0, 1])
        napari.run()"""

print("DONE!")