import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H2B+POLR3D'
file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 4
local_size = 200
rmax = 100


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


data = pd.DataFrame(columns=['nuclear', 'FOV', 'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear', 'mCherry_mean',
                             'mean_int_DNAFISH', 'n_ecDNA',
                             'mean_int_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'dis_to_hub_area',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized',
                             'radial_curve_DNAFISH_seg', 'nuclear_int',
                             'DNAFISH_int',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r',
                             'g', 'g0', 'g_value',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half'])

for fov in range(img_hoechst_stack.shape[0]):
    print("Analyzing %s, fov %s" % (sample, fov))
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry_bgc = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s%s/seg_tif/%s_%s_ecseg_n12.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

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
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_mCherry = img_mCherry_bgc.copy()
        local_mCherry = local_mCherry[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_nuclear_seg_mCherry = ima.img_local_seg(img_seg, position, nuclear_props[i].label)

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        local_IF_props = regionprops(label(local_nuclear_seg_mCherry), local_mCherry)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        perimeter_nuclear = local_nuclear_props[0].perimeter
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        mean_int_IF = local_IF_props[0].intensity_mean
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        circ_nuclear = (4 * math.pi * area_nuclear) / (perimeter_nuclear ** 2)

        # ecDNA measurements
        n_ecDNA = len(ecDNA_props)
        centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
        area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
        area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
        total_area_ecDNA = sum(area_ind_ecDNA)
        total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

        mean_int_ind_ecDNA = [ecDNA_props[i].intensity_mean for i in range(len(ecDNA_props))]
        total_int_ind_ecDNA = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA))
        total_int_ecDNA = sum(total_int_ind_ecDNA)
        mean_int_ecDNA = total_int_ecDNA / total_area_ecDNA if total_area_ecDNA != 0 else 0

        # percentage and cumulative curves
        percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True)) / total_area_ecDNA)
        percentage_area_ratio_ind_ecDNA = list(np.array(sorted(area_ratio_ind_ecDNA, reverse=True)) /
                                               total_area_ratio_ecDNA)
        percentage_total_int_ind_ecDNA = list(np.array(sorted(total_int_ind_ecDNA, reverse=True)) / total_int_ecDNA)

        cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
        cum_percentage_area_n_half = dat.find_pos(0.5, percentage_area_ind_ecDNA)
        cum_percentage_area_ratio_ind_ecDNA = dat.list_sum(percentage_area_ratio_ind_ecDNA)
        cum_percentage_area_ratio_n_half = dat.find_pos(0.5, percentage_area_ratio_ind_ecDNA)
        cum_percentage_total_int_ind_ecDNA = dat.list_sum(percentage_total_int_ind_ecDNA)
        cum_percentage_total_int_n_half = dat.find_pos(0.5, percentage_total_int_ind_ecDNA)

        cum_area_ind_ecDNA = dat.list_sum(area_ind_ecDNA)
        cum_area_n_half = dat.find_pos(cum_area_ind_ecDNA[-1] / 2, cum_area_ind_ecDNA)
        cum_area_ratio_ind_ecDNA = dat.list_sum(area_ratio_ind_ecDNA)
        cum_area_ratio_n_half = dat.find_pos(cum_area_ratio_ind_ecDNA[-1] / 2, cum_area_ratio_ind_ecDNA)
        cum_int_ind_ecDNA = dat.list_sum(total_int_ind_ecDNA)
        cum_int_n_half = dat.find_pos(cum_int_ind_ecDNA[-1] / 2, cum_int_ind_ecDNA)

        # distance from hub v2
        dis_to_hub_area_v2 = 0

        if n_ecDNA == 0:
            dis_to_hub_area_v2 = 0
        elif n_ecDNA > 1:
            ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'total_int': total_int_ind_ecDNA,
                                      'centroid': centroid_ind_ecDNA})

            ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False,
                                                               inplace=False,
                                                               ignore_index=True)
            dis = []
            for k in range(len(ind_ecDNA_sort_area)):
                dis_temp = 0
                target_ecDNA = ind_ecDNA_sort_area.iloc[[k]]
                rest_ecDNA = ind_ecDNA_sort_area.copy().drop(k)
                current_ecDNA = pd.concat([target_ecDNA, rest_ecDNA], axis=0).reset_index(drop=True)
                current_ecDNA['dis'] = \
                    [((current_ecDNA['centroid'][i][0] - current_ecDNA['centroid'][0][0]) ** 2 +
                      (current_ecDNA['centroid'][i][1] - current_ecDNA['centroid'][0][1]) ** 2) ** 0.5
                     for i in range(len(current_ecDNA))]

                for ind in range(n_ecDNA - 1):
                    dis_temp = dis_temp + current_ecDNA['area'][ind + 1] * current_ecDNA['dis'][
                        ind + 1] / total_area_ecDNA
                dis.append(dis_temp)
            ind_ecDNA_sort_area['dis'] = dis

            for ind in range(n_ecDNA):
                dis_to_hub_area_v2 = dis_to_hub_area_v2 + ind_ecDNA_sort_area['area'][ind] * \
                                     ind_ecDNA_sort_area['dis'][ind] / total_area_ecDNA

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                      local_DNAFISH, 0.1, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                      local_nuclear, 0.1, 1)
        radial_distribution_relative_r_DNAFISH_seg = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                      local_DNAFISH_seg, 0.1, 1)

        radial_distribution_relative_r_normalized = \
            list(np.array(radial_distribution_relative_r_DNAFISH) / np.array(
                radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_normalized = dat.nan_replace(radial_distribution_relative_r_normalized)

        nuclear_int = img_to_pixel_int(local_nuclear_seg, local_nuclear)
        DNAFISH_int = img_to_pixel_int(local_nuclear_seg, local_DNAFISH)
        DNAFISH_seg_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_r_to_center = img_to_pixel_int(local_nuclear_seg, local_centroid_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

        # auto-correlation
        _, r, g, dg = mat.auto_correlation(local_DNAFISH, local_nuclear_seg, rmax)
        g0 = g[0]
        g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

        data.loc[len(data.index)] = [i, fov, n_nuclear_convex_dilation,
                                     original_centroid_nuclear, area_nuclear, circ_nuclear, mean_int_nuclear, mean_int_IF,
                                     mean_int_DNAFISH, n_ecDNA,
                                     mean_int_ecDNA, total_area_ecDNA, total_area_ratio_ecDNA,
                                     dis_to_hub_area_v2,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                     radial_distribution_relative_r_normalized, radial_distribution_relative_r_DNAFISH_seg,
                                     nuclear_int, DNAFISH_int, DNAFISH_seg_label, int_r_to_edge, int_r_to_center, int_relative_r,
                                     g, g0, g_value,
                                     cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                     cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                     cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                     cum_area_ind_ecDNA, cum_area_n_half,
                                     cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half,
                                     cum_int_ind_ecDNA, cum_int_n_half]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])

data.to_csv('%s%s_n%s.txt' % (output_dir, sample, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")