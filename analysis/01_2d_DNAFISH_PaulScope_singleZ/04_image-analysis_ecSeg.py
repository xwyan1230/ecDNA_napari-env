from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import numpy as np
import shared.image as ima
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg
import tifffile as tif
import shared.math as mat
from skimage.filters import threshold_otsu
import shared.objects as obj
import matplotlib.pyplot as plt
import os
import napari

# INPUT PARAMETERS
# file info
# master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220816_Natasha_ColoDM_reimage/"
master_folder = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220825_POLR3D/20220825_POLR3Dtest/'
sample = 'Control_2'
master_path = '%s%s/' % (master_folder, sample)
version = 1
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100
rmax = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s_z.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

fovs = list(set(data_z['FOV'].tolist()))
total_fov = len(fovs)

data = pd.DataFrame(columns=['nuclear', 'FOV', 'z', 'z_min', 'z_max', 'z_ratio', 'label_nuclear', 'centroid_nuclear', 'limit',
                             'area_nuclear',
                             'mean_int_nuclear', 'total_int_nuclear',
                             'mean_int_DNAFISH', 'total_int_DNAFISH',
                             'n_ecDNA',
                             'mean_int_ecDNA', 'total_int_ecDNA',
                             'mean_int_ind_ecDNA', 'total_int_ind_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'max_area_ecDNA', 'max_area_ratio_ecDNA',
                             'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                             'dis_to_hub_area', 'dis_to_hub_int',
                             'g', 'g_value',
                             'radial_curve_DNAFISH', 'radial_curve_nuclear', 'radial_center', 'radial_edge',
                             'relative_r_area', 'relative_r_int',
                             'angle_curve_DNAFISH', 'angle_curve_nuclear', 'angle_value',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half'])

# IMAGING ANALYSIS
for f in range(total_fov):
    print("Analyzing %s, analyzing fov %s/%s" % (sample, f+1, total_fov))
    fov = fovs[f]
    data_z_fov = data_z[data_z['FOV'] == fov].copy().reset_index()

    # load images
    file_prefix = "%s_Position %s_RAW" % (sample, fov)
    im_z_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s_ch01.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_seg_convex = skio.imread("%s%s_seg.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_DNAFISH_seg_z = skio.imread("%s%s_ecSeg_z.tif" % (master_path, file_prefix), plugin="tifffile")

    for i in range(len(data_z_fov)):
        print("Analyzing %s, analyzing fov %s, nucleus %s/%s" % (sample, f+1, i+1, len(data_z_fov)))
        z_current = data_z_fov['z'][i]
        label_nuclear = data_z_fov['label_nuclear'][i]
        original_centroid_nuclear = data_z_fov['centroid_nuclear'][i]

        # get images for given z
        img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
        img_nuclear = im_z_stack_nuclear[z_current]
        img_DNAFISH = im_z_stack_DNAFISH[z_current]
        img_DNAFISH_seg = im_z_stack_DNAFISH_seg_z[z_current]

        # get local images
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_nuclear_centroid = local_nuclear_props[0].centroid
        local_DNAFISH_seg = img_DNAFISH_seg.copy()
        local_DNAFISH_seg = img_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg_convex == 0] = 0

        # background correction

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg_convex), local_DNAFISH)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        total_int_nuclear = area_nuclear * mean_int_nuclear

        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        total_int_DNAFISH = area_nuclear * mean_int_DNAFISH

        n_ecDNA = len(ecDNA_props)
        centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
        area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
        area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
        total_area_ecDNA = sum(area_ind_ecDNA)
        total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)
        max_area_ecDNA = max(area_ind_ecDNA + [0])
        max_area_ratio_ecDNA = max(area_ratio_ind_ecDNA + [0])

        mean_int_ind_ecDNA = [ecDNA_props[i].intensity_mean for i in range(len(ecDNA_props))]
        total_int_ind_ecDNA = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA))
        total_int_ecDNA = sum(total_int_ind_ecDNA)
        mean_int_ecDNA = total_int_ecDNA / total_area_ecDNA if total_area_ecDNA != 0 else 0

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
        cum_area_n_half = dat.find_pos(cum_area_ind_ecDNA[-1]/2, cum_area_ind_ecDNA)
        cum_area_ratio_ind_ecDNA = dat.list_sum(area_ratio_ind_ecDNA)
        cum_area_ratio_n_half = dat.find_pos(cum_area_ratio_ind_ecDNA[-1]/2, cum_area_ratio_ind_ecDNA)
        cum_int_ind_ecDNA = dat.list_sum(total_int_ind_ecDNA)
        cum_int_n_half = dat.find_pos(cum_int_ind_ecDNA[-1]/2, cum_int_ind_ecDNA)

        # distance from the hub
        dis_to_hub_area = 0
        dis_to_hub_int = 0
        if n_ecDNA == 0:
            dis_to_hub_area = -1
            dis_to_hub_int = -1
        elif n_ecDNA > 1:
            ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'total_int': total_int_ind_ecDNA, 'centroid': centroid_ind_ecDNA})

            ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
                                                               ignore_index=True)
            ind_ecDNA_sort_area['dis'] = \
                [((ind_ecDNA_sort_area['centroid'][i][0] - ind_ecDNA_sort_area['centroid'][0][0]) ** 2 +
                  (ind_ecDNA_sort_area['centroid'][i][1] - ind_ecDNA_sort_area['centroid'][0][1]) ** 2) ** 0.5
                 for i in range(len(ind_ecDNA_sort_area))]

            for ind in range(n_ecDNA - 1):
                dis_to_hub_area = dis_to_hub_area + ind_ecDNA_sort_area['area'][ind + 1] / total_area_ecDNA * \
                                  ind_ecDNA_sort_area['dis'][ind + 1]

            ind_ecDNA_sort_int = ind_ecDNA.copy().sort_values(by='total_int', axis=0, ascending=False,
                                                              inplace=False, ignore_index=True)
            ind_ecDNA_sort_int['dis'] = \
                [((ind_ecDNA_sort_int['centroid'][i][0] - ind_ecDNA_sort_int['centroid'][0][0]) ** 2 +
                  (ind_ecDNA_sort_int['centroid'][i][1] - ind_ecDNA_sort_int['centroid'][0][1]) ** 2) ** 0.5
                 for i in range(len(ind_ecDNA_sort_int))]

            for ind in range(n_ecDNA - 1):
                dis_to_hub_int = dis_to_hub_int + ind_ecDNA_sort_int['total_int'][ind + 1] / total_int_ecDNA * \
                                 ind_ecDNA_sort_int['dis'][ind + 1]

        # auto-correlation
        _, r, g, dg = mat.auto_correlation(local_DNAFISH, local_nuclear_seg_convex, rmax)
        g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg_convex, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg_convex,
                                                                  local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg_convex == 0] = 0
        local_edge_distance_map[local_nuclear_seg_convex == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                      local_DNAFISH, 0.01, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                      local_nuclear, 0.01, 1)

        radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
        radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

        radial_subtract = \
            np.array(radial_distribution_relative_r_DNAFISH_smooth) - np.array(radial_distribution_relative_r_nuclear_smooth)
        radial_subtract_center = np.mean(radial_subtract[0:40])
        radial_subtract_edge = np.mean(radial_subtract[40:80])

        # relative_r
        relative_r_area = 0
        relative_r_int = 0
        if n_ecDNA == 0:
            relative_r_area = -1
            relative_r_int = -1
        elif n_ecDNA >= 1:
            ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'total_int': total_int_ind_ecDNA,
                                      'centroid': centroid_ind_ecDNA})

            ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
                                                               ignore_index=True)
            ind_ecDNA_sort_area['dis'] = \
                [local_relative_r_map[int(ind_ecDNA_sort_area['centroid'][i][0])][
                     int(ind_ecDNA_sort_area['centroid'][i][1])]
                 for i in range(len(ind_ecDNA_sort_area))]

            for ind in range(n_ecDNA):
                relative_r_area = relative_r_area + ind_ecDNA_sort_area['area'][ind] / total_area_ecDNA * \
                                  ind_ecDNA_sort_area['dis'][ind]

            ind_ecDNA_sort_int = ind_ecDNA.copy().sort_values(by='total_int', axis=0, ascending=False,
                                                              inplace=False, ignore_index=True)
            ind_ecDNA_sort_int['dis'] = \
                [local_relative_r_map[int(ind_ecDNA_sort_int['centroid'][i][0])][
                     int(ind_ecDNA_sort_int['centroid'][i][1])]
                 for i in range(len(ind_ecDNA_sort_int))]

            for ind in range(n_ecDNA):
                relative_r_int = relative_r_int + ind_ecDNA_sort_int['total_int'][ind] / total_int_ecDNA * \
                                 ind_ecDNA_sort_int['dis'][ind]

        # angle distribution
        local_angle_map = ima.angle_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
        angle_distribution_DNAFISH = ima.radial_distribution_from_distance_map(local_nuclear_seg_convex,
                                                                               local_angle_map,
                                                                               local_DNAFISH, 1, 360)
        angle_distribution_nuclear = ima.radial_distribution_from_distance_map(local_nuclear_seg_convex,
                                                                               local_angle_map,
                                                                               local_nuclear, 1, 360)

        angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
        angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
        angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
            dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)
        angle_value = angle_distribution_DNAFISH_smooth_centered[179]

        data.loc[len(data.index)] = [i, data_z['FOV'][i], data_z['z'][i], data_z['z_min'][i], data_z['z_max'][i],
                                     data_z['z_ratio'][i], data_z['label_nuclear'][i], data_z['centroid_nuclear'][i],
                                     data_z['limit'][i],
                                     area_nuclear,
                                     mean_int_nuclear, total_int_nuclear,
                                     mean_int_DNAFISH, total_int_DNAFISH,
                                     n_ecDNA,
                                     mean_int_ecDNA, total_int_ecDNA,
                                     mean_int_ind_ecDNA, total_int_ind_ecDNA,
                                     total_area_ecDNA, total_area_ratio_ecDNA, max_area_ecDNA, max_area_ratio_ecDNA,
                                     area_ind_ecDNA, area_ratio_ind_ecDNA,
                                     dis_to_hub_area, dis_to_hub_int,
                                     g, g_value,
                                     radial_distribution_relative_r_DNAFISH_smooth,
                                     radial_distribution_relative_r_nuclear_smooth, radial_subtract_center,
                                     radial_subtract_edge,
                                     relative_r_area, relative_r_int,
                                     angle_distribution_DNAFISH_smooth_centered,
                                     angle_distribution_nuclear_smooth_centered, angle_value,
                                     cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                     cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                     cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                     cum_area_ind_ecDNA, cum_area_n_half, cum_area_ratio_ind_ecDNA,
                                     cum_area_ratio_n_half, cum_int_ind_ecDNA, cum_int_n_half]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])

data.to_csv('%s%s_v%s.txt' % (master_folder, sample, version), index=False, sep='\t')

print("DONE!")
