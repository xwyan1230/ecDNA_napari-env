import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.objects as obj
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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220822_sp8_POLR3Dko_JunDMandHSR_singleZ/"
sample = 'DM_Jun_POLR3Dko'
save_path = master_folder
start_fov = 0
total_fov = 30

# set parameters
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
n_nuclear_convex_dilation = 1
convex_conversion_threshold = 0.9
local_size = 200

# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, fov %s/%s" % (sample, fov+1, total_fov))
    if fov < 10:
        file_prefix = '20220822_sp8_POLR3Dko_JunDMandHSR_singleZ_%s_s0%s' % (sample, fov)
    else:
        file_prefix = '20220822_sp8_POLR3Dko_JunDMandHSR_singleZ_%s_s%s' % (sample, fov)

    img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_ch02.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_IF = skio.imread("%s%s/%s_ch01.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_nuclear_seg_convex = skio.imread("%s%s/%s_seg.tif" % (master_folder, sample, file_prefix), plugin="tifffile")

    # Perform nuclear segmentation
    """img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                       threshold=convex_conversion_threshold))
    img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
    tif.imwrite("%s%s/%s_seg.tif" % (master_folder, sample, file_prefix), img_nuclear_seg_convex)"""

    nuclear_props = regionprops(img_nuclear_seg_convex)

    img_DNAFISH_seg = np.zeros_like(img_DNAFISH)
    mean_int_DNAFISH = np.sum(img_DNAFISH * img_nuclear_seg_convex)/np.sum(img_nuclear_seg_convex)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov+1, i+1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, nuclear_props[i].label)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_IF = img_IF.copy()
        local_IF = local_IF[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid

        # ecDNA segmentation
        local_DNAFISH_singlet = local_DNAFISH.copy()
        local_DNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
        otsu_threshold_val_local_DNAFISH = threshold_otsu(local_DNAFISH_singlet)

        if otsu_threshold_val_local_DNAFISH == 0:
            print("skip due to no intensity in DNA FISH channel")
        else:
            threshold_min = otsu_threshold_val_local_DNAFISH + (img_DNAFISH.max() - otsu_threshold_val_local_DNAFISH) / 4
            threshold_min7 = otsu_threshold_val_local_DNAFISH + (img_DNAFISH.max() - otsu_threshold_val_local_DNAFISH) / 6
            print(threshold_min)

            FISH_seg_local = np.zeros_like(local_DNAFISH)

            for k in range(15):
                local = threshold_local(local_DNAFISH, 10 * k + 7)
                out = local_DNAFISH_singlet > local * 1.1
                out = binary_erosion(out)
                if k < 1:
                    out = binary_dilation(out)
                    out = obj.remove_small(out, 30)
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH)
                    print(len(out_props))
                    for j in range(len(out_props)):
                        temp = np.zeros_like(local_DNAFISH)
                        temp[out_label == out_props[j].label] = 1
                        temp_outer_edge = binary_dilation(temp, disk(50))
                        temp_outer_edge[temp == 1] = 0
                        mean_int_outer_edge = np.sum(local_DNAFISH * temp_outer_edge) / np.sum(temp_outer_edge)
                        if (out_props[j].intensity_mean / mean_int_outer_edge > 1.5) & (out_props[j].area > 5):
                            FISH_seg_local[out_label == out_props[j].label] = 1
                else:
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH_singlet)
                    for j in range(len(out_props)):
                        if out_props[j].area > 5:
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
            FISH_seg = obj.remove_small(FISH_seg, 12)
            """FISH_seg = binary_dilation(FISH_seg)
            FISH_seg = obj.remove_small(FISH_seg, 12)
            FISH_seg = binary_erosion(FISH_seg, disk(2))
            FISH_seg = binary_dilation(FISH_seg, disk(5))
            FISH_seg = binary_erosion(FISH_seg, disk(4))
            FISH_seg = obj.remove_small(FISH_seg, 12)"""

            viewer = napari.Viewer()
            # viewer.add_image(local_nuclear, blending='additive', colormap='blue')
            viewer.add_image(local_DNAFISH, blending='additive', colormap='green')
            # viewer.add_image(local_IF, blending='additive', colormap='magenta')
            viewer.add_image(FISH_seg_local, blending='additive')
            viewer.add_image(FISH_seg_watershed, blending='additive')
            viewer.add_image(FISH_seg, blending='additive')
            napari.run()

            img_DNAFISH_seg = ima.image_paste_to(img_DNAFISH_seg, FISH_seg,
                                                 [int(original_centroid_nuclear[0] - local_centroid[0]),
                                                  int(original_centroid_nuclear[1] - local_centroid[1])])

    tif.imwrite("%s%s/%s_ecSeg1.tif" % (master_folder, sample, file_prefix), img_DNAFISH_seg)

    viewer = napari.view_image(img_nuclear, blending='additive', colormap='blue')
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green')
    viewer.add_image(img_IF, blending='additive', colormap='magenta')
    plt.imsave('%s%s/%s_napari-img.tiff' % (save_path, sample, file_prefix), dis.blending(viewer))
    viewer.close()

    viewer1 = napari.view_image(img_nuclear_seg_convex, blending='additive', colormap='blue')
    viewer1.add_image(img_DNAFISH_seg, blending='additive', colormap='green')
    plt.imsave('%s%s/%s_napari-seg1.tiff' % (save_path, sample, file_prefix), dis.blending(viewer1))
    viewer1.close()

print("DONE!")

