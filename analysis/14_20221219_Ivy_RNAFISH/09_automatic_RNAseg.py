import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt
import shared.image as ima
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage import segmentation
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import os
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


sample = 'PC9'
# fov_lst = [fov_to_str(i) for i in np.arange(2, 37, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 7, 1)]  # Colo320DM
# fov_lst = [fov_to_str(i) for i in np.arange(1, 31, 1)]  # Colo320HSR
# fov_lst = ['MycI2_'+fov_to_str(i) for i in np.arange(1, 17, 1)]  # HCT116
# fov_lst = list(np.arange(1, 9, 1)) + list(np.arange(10, 13, 1)) + [fov_to_str(i) for i in np.arange(13, 16, 1)] + [fov_to_str(i) for i in np.arange(17, 26, 1)] + [fov_to_str(i) for i in np.arange(27, 36, 1)] # PC3
fov_lst = [fov_to_str(i) for i in np.arange(1, 2, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 8, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(9, 19, 1)]  # PC9

# set parameters
n_nuclear_convex_dilation = 0
local_size = 150

for fov in range(len(fov_lst)):
    print(fov)

    file_name = '%s_%s_Lng_SVCC_Processed001_RAW' % (sample, fov_lst[fov])
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_RNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_RNAFISH_seg = np.zeros_like(img_RNAFISH)
    img_nuclear_seg_convex = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov_lst[fov]), plugin="tifffile")
    img_nuclear_seg = np.zeros_like(img_RNAFISH)

    """viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive')
    shapes_add = viewer.add_shapes(name='nuclear seg add', ndim=2)
    shapes_seg = viewer.add_shapes(name='nuclear seg remove', ndim=2)
    napari.run()

    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_seg.data, 'remove', img_nuclear_seg_convex)
    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_nuclear_seg_convex)
    img_nuclear_seg_convex = obj.label_resort(img_nuclear_seg_convex)

    img_nuclear_seg = np.zeros_like(img_nuclear_seg_convex)
    # remove tiny borders
    nuclear_props_temp = regionprops(img_nuclear_seg_convex)
    for i in range(len(nuclear_props_temp)):
        if nuclear_props_temp[i].area > 500:
            img_nuclear_seg[img_nuclear_seg_convex == nuclear_props_temp[i].label] = nuclear_props_temp[i].label
    img_nuclear_seg = obj.label_resort(img_nuclear_seg)

    if not os.path.exists("%s%s/seg_tif_manual/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif_manual/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif_manual/%s_%s_seg.tif" % (output_dir, sample, sample, fov_lst[fov]), img_nuclear_seg_convex)

    img_nuclear_seg_convex = img_nuclear_seg.copy()"""

    # RNA FISH segmentation
    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
    nuclear_props = regionprops(img_nuclear_seg_convex)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, nuclear_props[i].label)
        local_RNAFISH = img_RNAFISH.copy()
        local_RNAFISH = local_RNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid

        local_RNAFISH_singlet = local_RNAFISH.copy()
        local_RNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
        otsu_threshold_val_local_RNAFISH = threshold_otsu(local_RNAFISH_singlet)

        if otsu_threshold_val_local_RNAFISH == 0:
            print("skip due to no intensity in RNA FISH channel")
        else:
            threshold_min = otsu_threshold_val_local_RNAFISH + (
                    img_RNAFISH.max() - otsu_threshold_val_local_RNAFISH) / 3
            threshold_min7 = otsu_threshold_val_local_RNAFISH + (
                    img_RNAFISH.max() - otsu_threshold_val_local_RNAFISH) / 6

            FISH_seg_local = np.zeros_like(local_RNAFISH)

            for k in range(15):
                local = threshold_local(local_RNAFISH, 10 * k + 7)
                out = local_RNAFISH_singlet > local
                out = binary_erosion(out)
                if k == 0:
                    out = binary_dilation(out)
                    out_label = label(out)
                    out_props = regionprops(out_label, local_RNAFISH_singlet)
                    for j in range(len(out_props)):
                        temp = np.zeros_like(local_RNAFISH)
                        temp[out_label == out_props[j].label] = 1
                        temp_outer_edge = binary_dilation(temp, disk(6))
                        temp_outer_edge[temp == 1] = 0
                        mean_int_outer_edge = np.sum(local_RNAFISH * temp_outer_edge) / np.sum(temp_outer_edge)
                        if (out_props[j].intensity_mean / mean_int_outer_edge > 1.2) & (out_props[j].area > 5) & \
                                (out_props[j].intensity_mean > threshold_min7):
                            FISH_seg_local[out_label == out_props[j].label] = 1
                else:
                    out_label = label(out)
                    out_props = regionprops(out_label, local_RNAFISH_singlet)
                    for j in range(len(out_props)):
                        if (out_props[j].intensity_mean > threshold_min) & (out_props[j].area > 5):
                            FISH_seg_local[out_label == out_props[j].label] = 1

            FISH_seg_watershed = np.zeros_like(local_RNAFISH)
            bg_val = otsu_threshold_val_local_RNAFISH * 3
            extreme_val = int(local_RNAFISH_singlet.max() * 2 / otsu_threshold_val_local_RNAFISH)
            maxima = extrema.h_maxima(local_RNAFISH, extreme_val)
            elevation_map = sobel(local_RNAFISH)
            markers = np.zeros_like(local_RNAFISH)
            markers[local_RNAFISH_singlet < bg_val] = 1
            markers[maxima == 1] = 2
            seg_wat = segmentation.watershed(elevation_map, markers)
            seg_wat_label = label(seg_wat)
            seg_wat_props = regionprops(seg_wat_label, local_RNAFISH_singlet)
            for j in range(len(seg_wat_props)):
                if (seg_wat_props[j].intensity_mean > threshold_min) & (seg_wat_props[j].area > 12):
                    FISH_seg_watershed[seg_wat_label == seg_wat_props[j].label] = 1

            FISH_seg = FISH_seg_watershed.copy()
            FISH_seg[FISH_seg_local == 1] = 1
            FISH_seg[local_nuclear_seg_convex == 0] = 0

            if FISH_seg.max() > 0:
                img_RNAFISH_seg = ima.image_paste_to(img_RNAFISH_seg, FISH_seg,
                                                 [int(original_centroid_nuclear[0] - local_centroid[0]),
                                                  int(original_centroid_nuclear[1] - local_centroid[1])])
                img_nuclear_seg[img_nuclear_seg_convex == nuclear_props[i].label] = nuclear_props[i].label
                img_nuclear_seg = obj.label_resort(img_nuclear_seg)

    if not os.path.exists("%s%s/seg_tif_filter/" % (output_dir, sample)):
        os.makedirs("%s%s/seg_tif_filter/" % (output_dir, sample))
    tif.imwrite("%s%s/seg_tif_filter/%s_%s_RNAseg.tif" % (output_dir, sample, sample, fov_lst[fov]), img_RNAFISH_seg)
    tif.imwrite("%s%s/seg_tif_filter/%s_%s_seg.tif" % (output_dir, sample, sample, fov_lst[fov]), img_nuclear_seg)

    if not os.path.exists("%s%s/color_img_filter/" % (output_dir, sample)):
        os.makedirs("%s%s/color_img_filter/" % (output_dir, sample))

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    plt.imsave('%s%s/color_img_filter/%s_%s_img.tiff' % (output_dir, sample, sample, fov_lst[fov]), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_RNAFISH_seg, blending='additive', colormap='red', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img_filter/%s_%s_seg.tiff' % (output_dir, sample, sample, fov_lst[fov]), dis.blending(viewer))
    viewer.close()

