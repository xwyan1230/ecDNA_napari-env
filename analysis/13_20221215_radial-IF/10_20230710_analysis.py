import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis, binary_erosion
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

sample = 'H3K27me2me3'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")
img_IF_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


# set parameters
n_nuclear_convex_dilation = 0
local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV', 'centroid_nuclear',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_IF',
                             'radial_curve_edge_nuclear', 'radial_curve_edge_DNAFISH', 'radial_curve_edge_IF',
                             'nuclear_int', 'DNAFISH_int', 'IF_int', 'DNAFISH_seg_label',
                             'int_r_to_edge', 'int_r_to_center', 'int_relative_r', 'mean_int_IF', 'mean_int_IF_ecDNA',
                             'mean_int_IF_nonecDNA'])

for fov in range(img_IF_stack.shape[0]):
    print(fov)
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_IF = img_IF_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]

    img_nuclear = np.concatenate([np.zeros(shape=[5, 3144]), img_nuclear], axis=0)[:3144, :3144]
    img_DNAFISH = np.concatenate([np.zeros(shape=[8, 3144]), img_DNAFISH], axis=0)[:3144, :3144]

    img_seg = skio.imread("%s%s/seg_tif_new/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s%s/seg_tif_new/%s_%s_ecseg1.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))
    else:
        img_nuclear_seg = img_seg

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
        local_DNAFISH_seg_erose = binary_erosion(local_DNAFISH_seg)
        local_nonDNAFISH_seg = local_nuclear_seg.copy()
        local_nonDNAFISH_seg[local_DNAFISH_seg_erose == 1] = 0

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_IF_props = regionprops(label(local_nuclear_seg), local_IF)
        local_IF_ecDNA_props = regionprops(label(local_DNAFISH_seg_erose), local_IF)
        local_IF_nonecDNA_props = regionprops(label(local_nonDNAFISH_seg), local_IF)

        mean_int_IF = local_IF_props[0].intensity_mean
        mean_int_IF_ecDNA_total_lst = [local_IF_ecDNA_props[i].intensity_mean * local_IF_ecDNA_props[i].area for i in range(len(local_IF_ecDNA_props))]
        mean_int_IF_ecDNA_area_lst = [local_IF_ecDNA_props[i].area for i in range(len(local_IF_ecDNA_props))]
        mean_int_IF_ecDNA = sum(mean_int_IF_ecDNA_total_lst)/(sum(mean_int_IF_ecDNA_area_lst)+0.001)
        mean_int_IF_nonecDNA_total_lst = [local_IF_nonecDNA_props[i].intensity_mean * local_IF_nonecDNA_props[i].area for i in
                                       range(len(local_IF_nonecDNA_props))]
        mean_int_IF_nonecDNA_area_lst = [local_IF_nonecDNA_props[i].area for i in range(len(local_IF_nonecDNA_props))]
        mean_int_IF_nonecDNA = sum(mean_int_IF_nonecDNA_total_lst) / (sum(mean_int_IF_nonecDNA_area_lst)+0.001)

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH, 0.025, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_nuclear, 0.025, 1)
        radial_distribution_relative_r_IF = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_IF, 0.025, 1)

        radial_distribution_edge_DNAFISH = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                     local_edge_distance_map,
                                                                                     local_DNAFISH, 1.5, 60)
        radial_distribution_edge_nuclear = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                     local_edge_distance_map,
                                                                                     local_nuclear, 1.5, 60)
        radial_distribution_edge_IF = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                     local_edge_distance_map,
                                                                                     local_IF, 1.5, 60)

        nuclear_int = img_to_pixel_int(local_nuclear_seg, local_nuclear)
        DNAFISH_int = img_to_pixel_int(local_nuclear_seg, local_DNAFISH)
        IF_int = img_to_pixel_int(local_nuclear_seg, local_IF)
        DNAFISH_seg_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_r_to_center = img_to_pixel_int(local_nuclear_seg, local_centroid_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

        data.loc[len(data.index)] = [i, fov, original_centroid_nuclear,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                     radial_distribution_relative_r_IF,
                                     radial_distribution_edge_nuclear, radial_distribution_edge_DNAFISH, radial_distribution_edge_IF,
                                     nuclear_int, DNAFISH_int, IF_int, DNAFISH_seg_label,
                                     int_r_to_edge, int_r_to_center, int_relative_r, mean_int_IF, mean_int_IF_ecDNA, mean_int_IF_nonecDNA]

data.to_csv('%s%s_radial_new_n%s.txt' % (output_dir, sample, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")