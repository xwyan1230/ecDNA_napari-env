import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import napari

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K4me3'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")
img_IF_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


# set parameters
n_nuclear_convex_dilation = 8
local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV', 'centroid_nuclear',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_IF',
                             'radial_curve_DNAFISH_normalized', 'radial_curve_IF_normalized',
                             'nuclear_int', 'DNAFISH_int', 'IF_int', 'ecDNA_label',
                             'int_r_to_edge', 'int_r_to_center', 'int_relative_r'])

for fov in range(img_IF_stack.shape[0]):
    print(fov)
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_IF = img_IF_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]

    img_nuclear = np.concatenate([np.zeros(shape=[3, 3144]), img_nuclear], axis=0)[:3144, :3144]
    img_DNAFISH = np.concatenate([np.zeros(shape=[6, 3144]), img_DNAFISH], axis=0)[:3144, :3144]

    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s%s/seg_tif/%s_%s_ecseg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    # measure
    # get local images
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
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH, 0.05, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_nuclear, 0.05, 1)
        radial_distribution_relative_r_IF = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_IF, 0.05, 1)

        radial_distribution_relative_r_DNAFISH_normalized = \
            list(np.array(radial_distribution_relative_r_DNAFISH) / np.array(radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_DNAFISH_normalized = dat.nan_replace(radial_distribution_relative_r_DNAFISH_normalized)

        radial_distribution_relative_r_IF_normalized = \
            list(np.array(radial_distribution_relative_r_IF) / np.array(radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_IF_normalized = dat.nan_replace(radial_distribution_relative_r_IF_normalized)

        nuclear_int = img_to_pixel_int(local_nuclear_seg, local_nuclear)
        DNAFISH_int = img_to_pixel_int(local_nuclear_seg, local_DNAFISH)
        IF_int = img_to_pixel_int(local_nuclear_seg, local_IF)
        ecDNA_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_r_to_center = img_to_pixel_int(local_nuclear_seg, local_centroid_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

        data.loc[len(data.index)] = [i, fov, original_centroid_nuclear,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                     radial_distribution_relative_r_IF,
                                     radial_distribution_relative_r_DNAFISH_normalized,
                                     radial_distribution_relative_r_IF_normalized,
                                     nuclear_int, DNAFISH_int, IF_int, ecDNA_label,
                                     int_r_to_edge, int_r_to_center, int_relative_r]

# data.to_csv('%s%s_radial_n%s.txt' % (output_dir, sample, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")