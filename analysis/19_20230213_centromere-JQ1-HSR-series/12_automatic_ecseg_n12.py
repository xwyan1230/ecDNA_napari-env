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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230206_analysis_centromere-JQ1-DM-HSR/"
data_dir1 = "%sdata/20230130_JQ1-series/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-JQ1_mix_mCh-Ctrl'
start_fov = 0
file_name = '%s_mChIF-DNAFISH_RAW' % sample
img_hoechst_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir1, sample, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir1, sample, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s/%s_ch02.tif" % (data_dir1, sample, file_name), plugin="tifffile")

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
n_nuclear_convex_dilation = 12
convex_conversion_threshold = 0.8
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

print(img_mCherry_stack.shape[0])

for f in range(img_mCherry_stack.shape[0]):
    fov = start_fov + f
    print("analyzing %s, fov: %s" % (sample, fov))
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    img_nuclear_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    # ecDNA segmentation
    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg_convex = dilation(img_nuclear_seg, disk(n_nuclear_convex_dilation))
    nuclear_props = regionprops(img_nuclear_seg_convex)

    img_DNAFISH_seg = np.zeros_like(img_DNAFISH)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, nuclear_props[i].label)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid

        local_DNAFISH_singlet = local_DNAFISH.copy()
        local_DNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
        otsu_threshold_val_local_DNAFISH = threshold_otsu(local_DNAFISH_singlet)

        if otsu_threshold_val_local_DNAFISH == 0:
            print("skip due to no intensity in DNA FISH channel")
        else:
            threshold_min = otsu_threshold_val_local_DNAFISH + (
                        img_DNAFISH.max() - otsu_threshold_val_local_DNAFISH) / 3
            threshold_min7 = otsu_threshold_val_local_DNAFISH + (
                        img_DNAFISH.max() - otsu_threshold_val_local_DNAFISH) / 6

            FISH_seg_local = np.zeros_like(local_DNAFISH)

            for k in range(15):
                local = threshold_local(local_DNAFISH, 10 * k + 7)
                out = local_DNAFISH_singlet > local
                out = binary_erosion(out)
                if k == 0:
                    out = binary_dilation(out)
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH_singlet)
                    for j in range(len(out_props)):
                        temp = np.zeros_like(local_DNAFISH)
                        temp[out_label == out_props[j].label] = 1
                        temp_outer_edge = binary_dilation(temp, disk(6))
                        temp_outer_edge[temp == 1] = 0
                        mean_int_outer_edge = np.sum(local_DNAFISH * temp_outer_edge) / np.sum(temp_outer_edge)
                        if (out_props[j].intensity_mean / mean_int_outer_edge > 1.2) & (out_props[j].area > 5) & \
                                (out_props[j].intensity_mean > threshold_min7):
                            FISH_seg_local[out_label == out_props[j].label] = 1
                else:
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH_singlet)
                    for j in range(len(out_props)):
                        if (out_props[j].intensity_mean > threshold_min) & (out_props[j].area > 5):
                            FISH_seg_local[out_label == out_props[j].label] = 1

            FISH_seg_watershed = np.zeros_like(local_DNAFISH)
            bg_val = otsu_threshold_val_local_DNAFISH * 3
            extreme_val = int(local_DNAFISH_singlet.max() * 2 / otsu_threshold_val_local_DNAFISH)
            maxima = extrema.h_maxima(local_DNAFISH, extreme_val)
            elevation_map = sobel(local_DNAFISH)
            markers = np.zeros_like(local_DNAFISH)
            markers[local_DNAFISH_singlet < bg_val] = 1
            markers[maxima == 1] = 2
            seg_wat = segmentation.watershed(elevation_map, markers)
            seg_wat_label = label(seg_wat)
            seg_wat_props = regionprops(seg_wat_label, local_DNAFISH_singlet)
            for j in range(len(seg_wat_props)):
                if (seg_wat_props[j].intensity_mean > threshold_min) & (seg_wat_props[j].area > 12):
                    FISH_seg_watershed[seg_wat_label == seg_wat_props[j].label] = 1

            FISH_seg = FISH_seg_watershed.copy()
            FISH_seg[FISH_seg_local == 1] = 1
            FISH_seg[local_nuclear_seg_convex == 0] = 0

            img_DNAFISH_seg = ima.image_paste_to(img_DNAFISH_seg, FISH_seg,
                                                 [int(original_centroid_nuclear[0] - local_centroid[0]),
                                                  int(original_centroid_nuclear[1] - local_centroid[1])])

    tif.imwrite("%s%s/seg_tif/%s_%s_ecseg_n%s.tif" % (output_dir, sample, sample, fov, n_nuclear_convex_dilation), img_DNAFISH_seg)

    if not os.path.exists("%s%s/color_img/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s_%s_seg_n%s.tiff' % (output_dir, sample, sample, fov, n_nuclear_convex_dilation), dis.blending(viewer))
    viewer.close()

print("DONE!")