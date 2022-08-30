import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.objects as obj
import imutils
import shared.image as ima
from skimage import segmentation
import numpy as np
import tifffile as tif
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import math
import napari
from scipy import ndimage
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/"
sample = 'DM'
sample_folder = '%s_singleZ_imageJ' % sample
save_path = master_folder

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
n_nuclear_convex_dilation = 1
convex_conversion_threshold = 0.9
local_size = 300

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

# load images
file_prefix = '%s_singleZ_25pos_RAW' % sample
img_nuclear_stack = skio.imread("%s%s/%s_ch00.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch02.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_IF_stack = skio.imread("%s%s/%s_ch01.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_nuclear_seg_stack = skio.imread("%s%s/%s_seg.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_DNAFISH_seg_stack = skio.imread("%s%s/%s_ecSeg.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")

total_fov = img_nuclear_stack.shape[0]
start_fov = 0

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, fov %s/%s" % (sample, fov+1, total_fov))
    avg_neighbour_overlay_DNAFISH = np.zeros(shape=[800, 800], dtype=np.uint16)
    avg_neighbour_overlay_nuclear = np.zeros(shape=[800, 800], dtype=np.uint16)
    avg_img_center = [400, 400]

    img_nuclear = img_nuclear_stack[fov]
    img_DNAFISH = img_DNAFISH_stack[fov]
    img_IF = img_IF_stack[fov]
    img_nuclear_seg = img_nuclear_seg_stack[fov]
    img_DNAFISH_seg = img_DNAFISH_seg_stack[fov]

    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_IF = img_IF.copy()
        local_IF = local_IF[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_DNAFISH_seg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_nuclear_props = regionprops(label(local_nuclear_seg))
        local_nuclear_centroid = local_nuclear_props[0].centroid

        # angle distribution
        local_angle_map = ima.angle_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        angle_distribution_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_angle_map, local_DNAFISH, 1, 360)
        max_angle = angle_distribution_DNAFISH.index(max(angle_distribution_DNAFISH))
        rotated_local_DNAFISH = imutils.rotate(local_DNAFISH, angle=max_angle)
        rotated_local_nuclear = imutils.rotate(local_nuclear, angle=max_angle)
        rotated_local_nuclear_seg = imutils.rotate(local_nuclear_seg, angle=max_angle)
        rotated_nuclear_props = regionprops(label(rotated_local_nuclear_seg))
        rotated_local_nuclear_centroid = rotated_nuclear_props[0].centroid

        direction = list((np.array(avg_img_center) - np.array(rotated_local_nuclear_centroid)).astype(int))

        mean_intensity_DNAFISH = np.mean(rotated_local_DNAFISH)
        avg_neighbour_overlay_DNAFISH = ima.sum_up_image(avg_neighbour_overlay_DNAFISH, rotated_local_DNAFISH,
                                                         direction, 5000.0 / mean_intensity_DNAFISH)

        mean_intensity_nuclear = np.mean(rotated_local_nuclear)
        avg_neighbour_overlay_nuclear = ima.sum_up_image(avg_neighbour_overlay_nuclear, rotated_local_nuclear,
                                                         direction, 5000.0 / mean_intensity_nuclear)

    viewer = napari.Viewer()
    viewer.add_image(avg_neighbour_overlay_nuclear, blending='additive', colormap='viridis')
    plt.imsave('%s%s/%s_napari_neighbour-overlay-nuclear_fov%s.tiff' % (save_path, sample_folder, sample, fov),
               dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(avg_neighbour_overlay_DNAFISH, blending='additive', colormap='viridis')
    plt.imsave('%s%s/%s_napari_neighbour-overlay-DNAFISH_fov%s.tiff' % (save_path, sample_folder, sample, fov),
               dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(avg_neighbour_overlay_nuclear, blending='additive', colormap='blue')
    viewer.add_image(avg_neighbour_overlay_DNAFISH, blending='additive', colormap='green')
    plt.imsave('%s%s/%s_napari_neighbour-overlay-img_fov%s.tiff' % (save_path, sample_folder, sample, fov),
               dis.blending(viewer))
    viewer.close()


print("DONE!")





