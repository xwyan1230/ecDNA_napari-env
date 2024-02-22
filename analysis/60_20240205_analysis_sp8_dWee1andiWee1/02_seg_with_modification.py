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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'DMSO_24hr'
file_name = '20240205_sp8_Wee1inhibitoranddegrader_acidFISH_DM_DMSO_pos9_1'
# file_name = '20240205_sp8_Wee1inhibitoranddegrader_acidFISH_dWee1_1uM_pos9_1'
total_fov = 9

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

for i in range(total_fov):
    img_nuclear = skio.imread("%s%s/%s_s%s_ch00.tif" % (data_dir, sample, file_name, i), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_s%s_ch01.tif" % (data_dir, sample, file_name, i), plugin="tifffile")

    # nuclear seg
    if os.path.exists("%s/seg_tif/%s/%s_%s_seg.tif" % (output_dir, sample, file_name, i)):
        img_nuclear_seg_convex = skio.imread("%s/seg_tif/%s/%s_%s_seg.tif" % (output_dir, sample, file_name, i), plugin="tifffile")
    else:
        img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)

        img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))

    # ecDNA segmentation
    """if os.path.exists("%s/seg_tif/%s/%s_%s_ecseg1.tif" % (output_dir, sample, file_name, i)):
        img_DNAFISH_seg1 = skio.imread("%s/seg_tif/%s/%s_%s_ecseg1.tif" % (output_dir, sample, file_name, i),
                                             plugin="tifffile")
    else:
        img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
        img_DNAFISH_seg1[img_DNAFISH > 10000] = 1
        img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 20)"""
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)
    img_DNAFISH_seg1[img_DNAFISH > 6000] = 1
    img_DNAFISH_seg1 = obj.remove_small(img_DNAFISH_seg1, 15)

    if not os.path.exists("%s/color_img/%s/" % (output_dir, sample)):
        os.makedirs("%s/color_img/%s/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 25000])
    plt.imsave('%s/color_img/%s/%s_%s_DNAFISH.tiff' % (output_dir, sample, file_name, i), dis.blending(viewer))
    viewer.close()

    img_nuclear_copy = img_nuclear.copy()
    img_nuclear_copy[img_nuclear < 5000] = 0
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_copy, blending='additive', colormap='blue', contrast_limits=[0, 45000])
    plt.imsave('%s/color_img/%s/%s_%s_nuclei.tiff' % (output_dir, sample, file_name, i), dis.blending(viewer))
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 25000])
    plt.imsave('%s/color_img/%s/%s_%s_img.tiff' % (output_dir, sample, file_name, i), dis.blending(viewer))
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='red', contrast_limits=[0, img_nuclear_seg_convex.max()])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    shapes_add = viewer.add_shapes(name='add', ndim=2)
    shapes_remove = viewer.add_shapes(name='remove', ndim=2)
    shapes_ecSeg_add = viewer.add_shapes(name='ecSeg add', ndim=2)
    shapes_ecSeg_remove = viewer.add_shapes(name='ecSeg remove', ndim=2)
    napari.run()

    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_remove.data, 'remove', img_nuclear_seg_convex)
    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_nuclear_seg_convex)
    img_nuclear_seg_convex = obj.label_resort(img_nuclear_seg_convex)

    img_DNAFISH_seg1 = ima.napari_add_or_remove(shapes_ecSeg_remove.data, 'remove', img_DNAFISH_seg1)
    img_DNAFISH_seg1 = ima.napari_add_or_remove(shapes_ecSeg_add.data, 'add_disk', img_DNAFISH_seg1)

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s/color_img/%s/%s_%s_seg.tiff' % (output_dir, sample, file_name, i), dis.blending(viewer))
    viewer.close()

    if not os.path.exists("%s/seg_tif/%s/" % (output_dir, sample)):
        os.makedirs("%s/seg_tif/%s/" % (output_dir, sample))
    tif.imwrite("%s/seg_tif/%s/%s_%s_seg.tif" % (output_dir, sample, file_name, i), img_nuclear_seg_convex)
    tif.imwrite("%s/seg_tif/%s/%s_%s_ecseg1.tif" % (output_dir, sample, file_name, i), img_DNAFISH_seg1)

print("DONE!")

