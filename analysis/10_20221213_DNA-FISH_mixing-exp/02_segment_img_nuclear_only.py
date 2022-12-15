import skimage.io as skio
import napari
import shared.segmentation as seg
import shared.objects as obj
import tifffile as tif
import matplotlib.pyplot as plt
import shared.display as dis
import math
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221017_mixing-test_mCherry-series_after-heating/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-H2B-mCherry'
folder_name = '20221012_mixing-test_mCherry-series'
file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
n_nuclear_convex_dilation = 0
convex_conversion_threshold = 0.9
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

print(img_mCherry_stack.shape[0])

for fov in range(img_mCherry_stack.shape[0]):
    print(fov)
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]

    # nuclear segmentation
    img_nuclear_seg = seg.nuclear_seg1(img_hoechst, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))
    if not os.path.exists("%s%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif/%s_%s_seg.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    if not os.path.exists("%s%s/color_img/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img/" % (output_dir, sample))
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
    plt.imsave('%s%s/color_img/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")