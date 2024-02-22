import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.objects as obj
from skimage.morphology import disk, dilation, medial_axis
import shared.dataframe as dat
import math
import shared.display as dis
from skimage.morphology import disk, dilation
import tifffile as tif
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import numpy as np
import pandas as pd
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
data_dir2 = "%sfigures/" % master_folder
data_dir3 = "%sdata/max_proj/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1
bg_DNAFISH_proj = 8005.90911194839

df = pd.read_csv('%s%s/%s_centroid_chosen.txt' % (data_dir2, sample, sample), na_values=['.'], sep='\t')

n_nuclear_convex_dilation = 4
local_size = 200


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


data = pd.DataFrame(columns=['FOV', 'z', 'label',
                             'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear',
                             'mean_int_DNAFISH', 'n_ecDNA',
                             'mean_int_ecDNA', 'area_ecDNA_individual', 'area_ratio_ecDNA_individual', 'mean_int_ecDNA_individual',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'dis_to_hub_area',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_laminB',
                             'radial_curve_edge_nuclear', 'radial_curve_edge_DNAFISH', 'radial_curve_edge_laminB',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half',
                             'n_ecDNA_lst', 'total_area_DNAFISH_lst', 'total_area_ratio_DNAFISH_lst', 'total_int_DNAFISH_proj',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r'])

for i in range(len(df)):
    fov = int(df['FOV'][i])
    z = int(df['z'][i])
    print('fov:%s z:%s' % (fov, z))
    if z < 10:
        file_name = '%s_%s_z0%s' % (prefix, fov, z)
    else:
        file_name = '%s_%s_z%s' % (prefix, fov, z)
    img_nuclear_bgc = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_laminB_bgc = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH_seg = skio.imread("%s%s/%s_ecseg1.tif" % (data_dir1, sample, file_name), plugin="tifffile")*np.ones_like(img_nuclear_bgc)
    img_seg = skio.imread("%s%s/%s_seg.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_DNAFISH_proj = skio.imread("%s%s/fov%s_DNAFISH_maxproj.tif" % (data_dir3, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    original_centroid_nuclear = [df['centroid_x'][i], df['centroid_y'][i]]
    position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
    local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, df['label'][i])
    local_nuclear = img_nuclear_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_laminB = img_laminB_bgc.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH_proj = img_DNAFISH_proj.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH_seg = img_DNAFISH_seg.copy()[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH_seg[local_nuclear_seg == 0] = 0

    # basic measurements
    local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
    local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
    local_DNAFISH_proj_props = regionprops(label(local_nuclear_seg), local_DNAFISH_proj)
    ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

    area_nuclear = local_nuclear_props[0].area
    perimeter_nuclear = local_nuclear_props[0].perimeter
    mean_int_nuclear = local_nuclear_props[0].intensity_mean
    mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
    circ_nuclear = (4 * math.pi * area_nuclear) / (perimeter_nuclear ** 2)

    # total ecDNA measurements
    total_int_DNAFISH_proj = (local_DNAFISH_proj_props[0].intensity_mean-bg_DNAFISH_proj) * area_nuclear

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

    radial_distribution_edge_DNAFISH = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                 local_edge_distance_map,
                                                                                 local_DNAFISH, 1.5, 60)
    radial_distribution_edge_laminB = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                 local_edge_distance_map,
                                                                                 local_laminB, 1.5, 60)
    radial_distribution_edge_nuclear = ima.radial_distribution_from_distance_map(local_nuclear_seg,
                                                                                 local_edge_distance_map,
                                                                                 local_nuclear, 1.5, 60)

    radial_distribution_relative_r_DNAFISH = \
        ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                  local_DNAFISH, 0.025, 1)
    radial_distribution_relative_r_laminB = \
        ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                  local_laminB, 0.025, 1)
    radial_distribution_relative_r_nuclear = \
        ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                  local_nuclear, 0.025, 1)



    n_ecDNA_lst = []
    total_area_DNAFISH_lst = []
    total_area_ratio_DNAFISH_lst = []
    for j in range(15):
        thresh = j * 5000
        local_DNAFISH_seg_temp = np.zeros_like(local_DNAFISH)
        local_DNAFISH_seg_temp[local_DNAFISH >= thresh] = 1
        local_DNAFISH_seg_temp = obj.remove_small(local_DNAFISH_seg_temp, 6)
        local_DNAFISH_seg_temp[local_nuclear_seg == 0] = 0
        ecDNA_props = regionprops(label(local_DNAFISH_seg_temp), local_DNAFISH)
        n_ecDNA_lst.append(len(ecDNA_props))

        local_DNAFISH_seg_temp = np.zeros_like(local_DNAFISH)
        local_DNAFISH_seg_temp[local_DNAFISH >= thresh] = 1
        local_DNAFISH_seg_temp[local_DNAFISH >= thresh + 5000] = 0
        local_DNAFISH_seg_temp[local_nuclear_seg == 0] = 0
        total_area_DNAFISH_lst.append(local_DNAFISH_seg_temp.sum())
        total_area_ratio_DNAFISH_lst.append(local_DNAFISH_seg_temp.sum() / area_nuclear)

    data.loc[len(data.index)] = [fov, z, df['label'][i],
                                 n_nuclear_convex_dilation,
                                 original_centroid_nuclear, area_nuclear, circ_nuclear, mean_int_nuclear,
                                 mean_int_DNAFISH, n_ecDNA,
                                 mean_int_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA, mean_int_ind_ecDNA,
                                 total_area_ecDNA, total_area_ratio_ecDNA,
                                 dis_to_hub_area_v2,
                                 radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                 radial_distribution_relative_r_laminB,
                                 radial_distribution_edge_nuclear, radial_distribution_edge_DNAFISH,
                                 radial_distribution_edge_laminB,
                                 cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                 cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                 cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                 cum_area_ind_ecDNA, cum_area_n_half,
                                 cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half,
                                 cum_int_ind_ecDNA, cum_int_n_half,
                                 n_ecDNA_lst, total_area_DNAFISH_lst, total_area_ratio_DNAFISH_lst, total_int_DNAFISH_proj,
                                 DNAFISH_seg_label, int_r_to_edge, int_relative_r]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])

if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
data.to_csv('%s%s/%s_n%s_full.txt' % (output_dir, sample, sample, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")