import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240205_analysis_sp8_dWee1andiWee1/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'iWee1_62point5nM_24hr'
sets = 1
file_name = '20240315_sp8_dWee1andiWee1_series_iWee1_62point5nM_24hr_9pos_%s' % sets
# file_name = '20240307_sp8_iWee1anddWee1_acidFISH_24hr_iWee1_125nM_24hr_9pos_%s' % sets
# file_name = '20240315_sp8_dWee1andiWee1_series_dWee1_500nM_24hr_9pos_%s' % sets
# file_name = '20240307_sp8_iWee1anddWee1_acidFISH_24hr_iWee1_1uM_24hr_9pos_%s' % sets
# file_name = '20240315_sp8_dWee1andiWee1_series_dWee1_2uM_24hr_9pos_%s' % sets
# file_name = '20240307_sp8_iWee1anddWee1_acidFISH_24hr_dWee1_125nM_24hr_4pos_%s' % sets  # set 2 fov 1
# file_name = '20240202_Colo320DM_dWee1andiWee1_acidFISH_Colo320DM_iWee1_250nM_pos9_%s' % sets
# file_name = '20240202_Colo320DM_dWee1andiWee1_acidFISH_Colo320DM_dWee1_250nM_pos9_%s' % sets
# file_name = '20240205_sp8_Wee1inhibitoranddegrader_acidFISH_iWee1_500nM_pos9_%s' % sets
# file_name = '20240205_sp8_Wee1inhibitoranddegrader_acidFISH_DM_DMSO_pos9_%s' % sets
# file_name = '20240205_sp8_Wee1inhibitoranddegrader_acidFISH_dWee1_1uM_pos9_%s' % sets
total_fov = 9
n_nuclear_convex_dilation = 0
local_size = 200


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


data = pd.DataFrame(columns=['sample', 'file_name', 'FOV', 'nuclear', 'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_DNAFISH', 'n_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'dis_to_hub_area',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r'])

for fov in range(total_fov):
    print("Analyzing fov %s" % (fov+1))
    img_nuclear = skio.imread("%s%s/%s_s%s_ch00.tif" % (data_dir, sample, file_name, fov), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_s%s_ch01.tif" % (data_dir, sample, file_name, fov), plugin="tifffile")
    img_seg = skio.imread("%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir1, sample, file_name, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s/seg_tif/%s/%s_%s_ecseg1.tif" % (data_dir1, sample, file_name, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))
    else:
        img_nuclear_seg = img_seg

    # measure
    # get local images
    nuclear_props = regionprops(img_nuclear_seg)
    img_nuclear_bgc = img_nuclear
    img_DNAFISH_bgc = img_DNAFISH

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov+1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        perimeter_nuclear = local_nuclear_props[0].perimeter
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
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean

        # radial measurements
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)
        DNAFISH_seg_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

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

        data.loc[len(data.index)] = [sample, file_name, fov, i, n_nuclear_convex_dilation,
                                     original_centroid_nuclear, area_nuclear, circ_nuclear, mean_int_DNAFISH, n_ecDNA,
                                     total_area_ecDNA, total_area_ratio_ecDNA, dis_to_hub_area_v2,
                                     DNAFISH_seg_label, int_r_to_edge, int_relative_r]

if not os.path.exists("%s/txt/" % output_dir):
    os.makedirs("%s/txt/" % output_dir)
data.to_csv('%s/txt/%s_%s_n%s.txt' % (output_dir, sample, sets, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")