import skimage.io as skio
import napari
import shared.segmentation as seg
import shared.objects as obj
import tifffile as tif
import matplotlib.pyplot as plt
import shared.display as dis
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage import segmentation
import math
import shared.image as ima
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K27me2me3'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir, file_name), plugin="tifffile")

img_IF_stack = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 0.9]  # used to filter nucleus
n_nuclear_convex_dilation = 2
convex_conversion_threshold = 0.8
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

for fov in range(img_IF_stack.shape[0]):
    print(fov)
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_IF = img_IF_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]

    img_nuclear = np.concatenate([np.zeros(shape=[3, 3144]), img_nuclear], axis=0)[:3144, :3144]
    img_DNAFISH = np.concatenate([np.zeros(shape=[6, 3144]), img_DNAFISH], axis=0)[:3144, :3144]

    # nuclear segmentation
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    napari.run()"""

    if not os.path.exists("%s%s/seg_tif_new/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif_new/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif_new/%s_%s_seg.tif" % (output_dir, sample, sample, fov), img_nuclear_seg_convex)

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 22000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 10)
    tif.imwrite("%s%s/seg_tif_new/%s_%s_ecseg1.tif" % (output_dir, sample, sample, fov), img_DNAFISH_seg1)

    if not os.path.exists("%s%s/color_img_new/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img_new/" % (output_dir, sample))
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    # viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%s%s/color_img_new/%s_%s_img.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_IF, blending='additive', colormap='red', contrast_limits=[0, img_IF.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%s%s/color_img_new/%s_%s_img_with_IF.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img_new/%s_%s_seg.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")