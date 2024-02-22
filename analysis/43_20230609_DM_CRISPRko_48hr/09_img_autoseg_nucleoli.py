import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import math
from skimage.measure import label, regionprops
from skimage.filters import try_all_threshold
from skimage.morphology import binary_dilation, binary_erosion, dilation, erosion, disk
import shared.objects as obj
import tifffile as tif
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'C'
sample = 'C2'
total_fov = 16
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
convex_conversion_threshold = 0.85
local_size = 200
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi

for fov in range(total_fov):
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
        # filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear_bgc = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_nuclear_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_seg.tif" % (data_dir, sample, sample, fov),
                          plugin="tifffile")

    nuclear_props = regionprops(img_nuclear_seg)
    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        label_nuclear = nuclear_props[i].label
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear_seg = binary_erosion(local_nuclear_seg, disk(3))
        local_nuclear = img_nuclear_bgc.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()[position[0]:position[1], position[2]:position[3]]

        local_nucleoli = np.zeros_like(local_nuclear)
        local_nucleoli[local_nuclear < 7000] = 1
        local_nucleoli[local_nuclear_seg == 0] = 0
        local_nucleoli_final = obj.remove_small(label(local_nucleoli), 10)
        local_nucleoli_final = dilation(local_nucleoli_final)
        local_nucleoli_final = ndimage.binary_fill_holes(local_nucleoli_final)
        local_nucleoli_final = obj.label_resort(obj.remove_small(local_nucleoli_final, 50))
        local_nucleoli_final = dilation(local_nucleoli_final)

        local_nucleoli_convex = obj.label_resort(seg.obj_to_convex_filter(label(local_nucleoli_final), threshold=0.5))

        viewer = napari.Viewer()
        viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, local_nuclear.max()])
        viewer.add_image(local_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, local_DNAFISH.max()])
        viewer.add_image(local_nucleoli_convex, blending='additive', contrast_limits=[0, 1])
        # plt.imsave('%sseg/%s/color_img/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
        napari.run()

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%sseg/%s/color_img/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()"""

print("DONE!")