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
start_fov = 0

for f in range(total_fov):
    fov = f+start_fov
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_nuclear_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_seg.tif" % (data_dir, sample, sample, fov), plugin="tifffile")
    img_ecDNA_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_ecseg1.tif" % (data_dir, sample, sample, fov), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_nuclear_seg, blending='additive', contrast_limits=[0, 1])
    shapes_add = viewer.add_shapes(name='add', ndim=2)
    shapes_remove = viewer.add_shapes(name='remove', ndim=2)
    napari.run()

    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_remove.data, 'remove', img_nuclear_seg)
    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_nuclear_seg_convex)
    img_nuclear_seg_convex = obj.label_resort(img_nuclear_seg_convex)

    if not os.path.exists("%sseg1/%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%sseg1/%s/seg_tif/" % (output_dir, sample))
    tif.imwrite("%sseg1/%s/seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    if not os.path.exists("%sseg1/%s/color_img/" % (output_dir, sample)):
        os.makedirs("%sseg1/%s/color_img/" % (output_dir, sample))
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%sseg1/%s/color_img/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_ecDNA_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%sseg1/%s/color_img/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    napari.run()

print("DONE!")