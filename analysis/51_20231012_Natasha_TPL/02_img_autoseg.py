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
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
from skimage.filters import threshold_otsu, threshold_local, threshold_yen, sobel
import tifffile as tif
import shared.segmentation as seg
import os
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_Natasha_colcemid/TPL/mCh-DMSO_GFP-TPL/"
output_dir = "%s04_figures/" % master_folder

sample = 'mCh-DMSO_GFP-TPL'
total_fov = 49
start_fov = 0
convex_conversion_threshold = 0.85
local_size = 200
frame = 4

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    if fov < 10:
        filename = '230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1_s0%s' % fov
    else:
        filename = '230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1_s%s' % fov
    img_nuclear = skio.imread("%s/03_DNAFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1/%s_ch00.tif" % (master_folder, filename), plugin="tifffile")
    img_DNAFISH = skio.imread("%s/03_DNAFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1/%s_ch01.tif" % (master_folder, filename), plugin="tifffile")

    # nuclear seg
    img_nuclear_seg = seg.nuclear_seg_global(img_nuclear, min_size=5000, max_size=30000)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg, threshold=convex_conversion_threshold))

    if not os.path.exists("%s/04_seg/%s/seg_tif/" % (master_folder, frame)):
        os.makedirs("%s/04_seg/%s/seg_tif/" % (master_folder, frame))
    tif.imwrite("%s/04_seg/%s/seg_tif/%s_%s_seg.tif" % (master_folder, frame, sample, fov), img_nuclear_seg_convex)

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 12000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 6)
    tif.imwrite("%s/04_seg/%s/seg_tif/%s_%s_ecseg1.tif" % (master_folder, frame, sample, fov), img_DNAFISH_seg1)

    if not os.path.exists("%s/04_seg/%s/color_img/" % (master_folder, frame)):
        os.makedirs("%s/04_seg/%s/color_img/" % (master_folder, frame))

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%s/04_seg/%s/color_img/%s_%s_img.tiff' % (master_folder, frame, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s/04_seg/%s/color_img/%s_%s_seg.tiff' % (master_folder, frame, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")