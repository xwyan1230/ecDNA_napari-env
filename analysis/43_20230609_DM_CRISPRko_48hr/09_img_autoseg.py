import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import math
import shared.objects as obj
import tifffile as tif
import shared.segmentation as seg
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'E2'
sample = 'E9'
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
        # filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
        filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        # filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
        filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")

    # nuclear seg
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    if not os.path.exists("%sseg/%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%sseg/%s/seg_tif/" % (output_dir, sample))
    tif.imwrite("%sseg/%s/seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 12000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 20)
    tif.imwrite("%sseg/%s/seg_tif/%s_%s_ecseg1.tif" % (output_dir, sample, sample, fov), img_DNAFISH_seg1)

    if not os.path.exists("%sseg/%s/color_img/" % (output_dir, sample)):
        os.makedirs("%sseg/%s/color_img/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%sseg/%s/color_img/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%sseg/%s/color_img/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")