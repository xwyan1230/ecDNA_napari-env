import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import math
import shared.objects as obj
import tifffile as tif
import imutils
import pandas as pd
import shared.segmentation as seg
import os

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8-1'
check_seg = 'NO'  # only accepts YES or NO

# NO NEED TO CHANGE
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
# treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s' % (sample, treatment, sample)
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s' % sample

total_fov = 16
start_fov = 0
pixel_size = 58.7  # 180  # nm
# 766nm for Keyence 10x
# 360nm for 512x512 at 63x
# 58.6nm for 3144x3144 at 63x (0.0765)
# 180nm for 1024x1024 at 63x (0.2349)
# 720nm for 256x256 at 63x (0.9399)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
convex_conversion_threshold = 0.8  # 0.8
maxima_threshold = 10  # 3
DNAFISH_threshold = 8000  # 8000
DNAFISH_remove_small_threshold = 15  # 3
local_factor_nuclear = 99  # 299  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi

if os.path.exists("%s/segmentation_parameters.txt" % output_dir):
    param = pd.read_csv('%s/segmentation_parameters.txt' % output_dir, na_values=['.'], sep='\t')
else:
    param = pd.DataFrame(columns=['sample', 'min_size_nuclear', 'max_size_nuclear', 'maxima_threshold', 'local_factor', 'convex_conversion_threshold',
                             'DNAFISH_threshold', 'DNAFISH_remove_small_threshold'])

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    img_nuclear_ori = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov + 1), plugin="tifffile")
    img_DNAFISH_ori = skio.imread("%s/DNAFISH/%s/%s_%s_ch00.tif" % (data_dir, sample, filename, fov + 1), plugin="tifffile")
    img_nuclear = cv2.flip(imutils.rotate(img_nuclear_ori, angle=-90), 0)
    img_DNAFISH = cv2.flip(imutils.rotate(img_DNAFISH_ori, angle=-90), 0)

    # nuclear seg
    if pixel_size == 180:
        img_nuclear_seg = seg.nuclear_seg_screen(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                         max_size=max_size_nuclear, maxima_threshold=maxima_threshold)
    else:
        img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                             max_size=max_size_nuclear, maxima_threshold=maxima_threshold)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    if not os.path.exists("%s/%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%s/%s/seg_tif/" % (output_dir, sample))
    tif.imwrite("%s/%s/seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov+1), img_nuclear_seg_convex)

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > DNAFISH_threshold] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, DNAFISH_remove_small_threshold)
    tif.imwrite("%s/%s/seg_tif/%s_%s_ecseg.tif" % (output_dir, sample, sample, fov+1), img_DNAFISH_seg1)

    # check segmentation
    if check_seg == 'YES':
        viewer = napari.Viewer()
        viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_seg_convex.max()])
        viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
        napari.run()

    if not os.path.exists("%s/%s/color_img/" % (output_dir, sample)):
        os.makedirs("%s/%s/color_img/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
    plt.imsave('%s/%s/color_img/%s_%s_img.tiff' % (output_dir, sample, sample, fov+1), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s/%s/color_img/%s_%s_seg.tiff' % (output_dir, sample, sample, fov+1), dis.blending(viewer))
    viewer.close()

if sample in param['sample'].tolist():
    sample_index = param[param['sample'] == sample].index[0]
    param.loc[sample_index] = [sample, min_size_nuclear, max_size_nuclear, maxima_threshold, local_factor_nuclear,
                                   convex_conversion_threshold, DNAFISH_threshold, DNAFISH_remove_small_threshold]
else:
    param.loc[len(param.index)] = [sample, min_size_nuclear, max_size_nuclear, maxima_threshold, local_factor_nuclear,
                                   convex_conversion_threshold, DNAFISH_threshold, DNAFISH_remove_small_threshold]
param.to_csv('%s/segmentation_parameters.txt' % output_dir, index=False, sep='\t')

print("DONE!")