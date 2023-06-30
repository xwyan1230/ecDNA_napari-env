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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
batch = 2
img_hoechst_stack = skio.imread("%s%s/%s_%s_%s_50pos_RAW_ch00.tif" % (data_dir, sample, exp, sample, batch), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_%s_%s_50pos_RAW_ch01.tif" % (data_dir, sample, exp, sample, batch), plugin="tifffile")

start_fov = 16

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
convex_conversion_threshold = 0.85
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000 / (pixel_size * 2)) ** 2 * math.pi

for f in range(img_hoechst_stack.shape[0]):
# for f in range(46):
    fov = start_fov + f
    print("analyzing %s, fov: %s" % (sample, fov))
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]

    # nuclear seg
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)

    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    if not os.path.exists("%s%s/seg_tif/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif/%s_%s_%s_seg.tif" % (output_dir, sample, sample, batch, fov), img_nuclear_seg_convex)

    img_nuclear_seg_convex_original = img_nuclear_seg_convex

    # ecDNA segmentation
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 12000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 20)

    tif.imwrite("%s%s/seg_tif/%s_%s_%s_ecseg1.tif" % (output_dir, sample, sample, batch, fov),img_DNAFISH_seg1)

    if not os.path.exists("%s%s/color_img/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%s%s/color_img/%s_%s_%s_img.tiff' % (output_dir, sample, sample, batch, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex_original, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_%s_%s_seg.tiff' % (output_dir, sample, sample, batch, fov), dis.blending(viewer))
    viewer.close()

print("DONE!")

