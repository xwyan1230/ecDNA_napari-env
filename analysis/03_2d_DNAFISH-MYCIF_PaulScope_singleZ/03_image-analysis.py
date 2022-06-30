from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import numpy as np
import shared.image as img
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg
import tifffile as tif
import matplotlib.pyplot as plt
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
sample_row = 'E'
sample_num = '10'
sample = 'E10'
raw_folder = '01_raw'
seg_folder = '02_seg'
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s/%s/%s_z.txt' % (master_folder, sample_row, sample_num, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

data = pd.DataFrame(columns=['nuclear',
                             'FOV',
                             'z',
                             'z_min',
                             'z_max',
                             'z_ratio',
                             'label_nuclear',
                             'centroid_nuclear',
                             'limit',
                             'area_nuclear',
                             'bg_int',
                             'mean_int_DNAFISH',
                             'mean_int_DNAFISH_norm',
                             'mean_int_nuclear',
                             'total_int_DNAFISH',
                             'total_int_DNAFISH_norm',
                             'total_int_nuclear',
                             'radial_curve_DNAFISH',
                             'radial_curve_DNAFISH_bg_correct',
                             'radial_curve_nuclear',
                             'radial_center',
                             'radial_edge',
                             'angle_curve_DNAFISH',
                             'angle_curve_DNAFISH_bg_correct',
                             'angle_curve_nuclear',
                             'total_area_ecDNA',
                             'area_ratio_ecDNA',
                             'total_int_ecDNA',
                             'total_int_ecDNA_norm',
                             'mean_int_ecDNA',
                             'mean_int_ecDNA_norm',
                             'area_individual_ecDNA',
                             'area_ratio_individual_ecDNA',
                             'mean_int_individual_ecDNA',
                             'mean_int_individual_ecDNA_norm',
                             'total_int_individual_ecDNA',
                             'total_int_individual_ecDNA_norm',
                             'percentage_area_curve_ecDNA',
                             'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA',
                             'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA',
                             'percentage_int_n_half',
                             'percentage_int_curve_ecDNA_norm',
                             'percentage_int_norm_n_half',
                             'cum_area_ind_ecDNA',
                             'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA',
                             'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA',
                             'cum_int_n_half',
                             'cum_int_ind_ecDNA_norm',
                             'cum_int_norm_n_half'])

# IMAGING ANALYSIS
for i in range(len(data_z)):
    print("Analyzing nucleus %s/%s" % (i+1, len(data_z)))
    fov = data_z['FOV'][i]
    z_current = data_z['z'][i]
    label_nuclear = data_z['label_nuclear'][i]
    original_centroid_nuclear = data_z['centroid_nuclear'][i]

    # load images
    im_z_stack_nuclear = skio.imread("%s%s/%s/%s/R%s_RAW_ch00.tif" %
                                     (master_folder, sample_row, sample_num, raw_folder, fov + 1), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s/R%s_RAW_ch01.tif" %
                                     (master_folder, sample_row, sample_num, raw_folder, fov + 1), plugin="tifffile")

    total_z = im_z_stack_nuclear.shape[0]
    im_z_stack_seg_convex = skio.imread("%s%s/%s/%s/R%s_seg.tif" %
                                        (master_folder, sample_row, sample_num, seg_folder, fov + 1), plugin="tifffile")

    # get images for given z
    img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
    img_nuclear = im_z_stack_nuclear[z_current]
    img_DNAFISH = im_z_stack_DNAFISH[z_current]

    # get local images
    position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
    local_nuclear = img_nuclear.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]

    # ecDNA segmentation
    int_thresh = data_z['limit'][i]*0.9
    k_dots = 5000
    vector = []
    vector_cum_weight = []
    weight = 0
    for m in range(len(local_nuclear_seg_convex)):
        for n in range(len(local_nuclear_seg_convex[0])):
            if local_nuclear_seg_convex[m][n] == 1:
                vector.append([m, n])
                if local_DNAFISH[m][n] > int_thresh:
                    weight = weight + local_DNAFISH[m][n] - int_thresh
                vector_cum_weight.append(weight)
    if weight != 0:
        random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=k_dots)
        img_dot = np.zeros_like(local_nuclear_seg_convex)
        for m in random_dot:
            img_dot[m[0]][m[1]] = img_dot[m[0]][m[1]] + 1
        img_dot_remove_bg = dilation(img_dot)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = dilation(img_dot_remove_bg)
        img_dot_seg = img_dot_remove_bg.copy()
        img_dot_seg[img_dot_remove_bg > 0] = 1
        DNAFISH_seg, _ = seg.find_organelle(local_DNAFISH, 'na', extreme_val=500, bg_val=data_z['limit'][i]*0.8,
                                            min_size=0, max_size=10000)
        DNAFISH_seg[local_nuclear_seg_convex == 0] = 0
        ecDNA_seg = img_dot_seg.copy()
        ecDNA_seg[DNAFISH_seg == 1] = 1
        bg_seg = local_nuclear_seg_convex.copy()
        bg_seg[ecDNA_seg == 1] = 0

        # background correction
        bg_local_DNAFISH = local_DNAFISH.copy()
        bg_local_DNAFISH[bg_seg == 0] = 0
        mean_int_bg = np.sum(bg_local_DNAFISH) / np.sum(bg_seg)
        local_DNAFISH_bg_correct = local_DNAFISH.copy()
        local_DNAFISH_bg_correct = local_DNAFISH_bg_correct.astype(float) - mean_int_bg
        local_DNAFISH_bg_correct[local_DNAFISH_bg_correct < 0] = 0

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg_convex), local_DNAFISH)
        ecDNA_props = regionprops(label(ecDNA_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        total_int_nuclear = area_nuclear * mean_int_nuclear

        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        total_int_DNAFISH = area_nuclear * mean_int_DNAFISH
        mean_int_DNAFISH_norm = mean_int_DNAFISH - mean_int_bg
        total_int_DNAFISH_norm = mean_int_DNAFISH_norm * area_nuclear

        n_ecDNA = len(ecDNA_props)
        area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
        area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
        total_area_ecDNA = sum(area_ind_ecDNA)
        total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

        mean_int_ind_ecDNA = [ecDNA_props[i].intensity_mean for i in range(len(ecDNA_props))]
        total_int_ind_ecDNA = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA))
        total_int_ecDNA = sum(total_int_ind_ecDNA)
        mean_int_ecDNA = total_int_ecDNA / total_area_ecDNA

        mean_int_ind_ecDNA_norm = list(np.array(mean_int_ind_ecDNA) - mean_int_bg)
        total_int_ind_ecDNA_norm = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA_norm))
        total_int_ecDNA_norm = sum(total_int_ind_ecDNA_norm)
        mean_int_ecDNA_norm = total_int_ecDNA_norm / total_area_ecDNA

        percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True)) / total_area_ecDNA)
        percentage_area_ratio_ind_ecDNA = list(np.array(sorted(area_ratio_ind_ecDNA, reverse=True)) /
                                               total_area_ratio_ecDNA)
        percentage_total_int_ind_ecDNA = list(np.array(sorted(total_int_ind_ecDNA, reverse=True)) / total_int_ecDNA)
        percentage_total_int_ind_ecDNA_norm = \
            list(np.array(sorted(total_int_ind_ecDNA_norm, reverse=True)) / total_int_ecDNA_norm)

        cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
        cum_percentage_area_n_half = dat.find_pos(0.5, percentage_area_ind_ecDNA)
        cum_percentage_area_ratio_ind_ecDNA = dat.list_sum(percentage_area_ratio_ind_ecDNA)
        cum_percentage_area_ratio_n_half = dat.find_pos(0.5, percentage_area_ratio_ind_ecDNA)
        cum_percentage_total_int_ind_ecDNA = dat.list_sum(percentage_total_int_ind_ecDNA)
        cum_percentage_total_int_n_half = dat.find_pos(0.5, percentage_total_int_ind_ecDNA)
        cum_percentage_total_int_ind_ecDNA_norm = dat.list_sum(percentage_total_int_ind_ecDNA_norm)
        cum_percentage_total_int_norm_n_half = dat.find_pos(0.5, percentage_total_int_ind_ecDNA_norm)

        cum_area_ind_ecDNA = dat.list_sum(area_ind_ecDNA)
        cum_area_n_half = dat.find_pos(cum_area_ind_ecDNA[-1]/2, cum_area_ind_ecDNA)
        cum_area_ratio_ind_ecDNA = dat.list_sum(area_ratio_ind_ecDNA)
        cum_area_ratio_n_half = dat.find_pos(cum_area_ratio_ind_ecDNA[-1]/2, cum_area_ratio_ind_ecDNA)
        cum_int_ind_ecDNA = dat.list_sum(total_int_ind_ecDNA)
        cum_int_n_half = dat.find_pos(cum_int_ind_ecDNA[-1]/2, cum_int_ind_ecDNA)
        cum_int_ind_ecDNA_norm = dat.list_sum(total_int_ind_ecDNA_norm)
        cum_int_norm_n_half = dat.find_pos(cum_int_ind_ecDNA_norm[-1]/2, cum_int_ind_ecDNA_norm)

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg_convex, return_distance=True)
        local_centroid_distance_map = img.distance_map_from_point(local_nuclear_seg_convex,
                                                                  local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg_convex == 0] = 0
        local_edge_distance_map[local_nuclear_seg_convex == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                      local_DNAFISH, 0.01, 1)
        radial_distribution_relative_r_DNAFISH_bg_correct = \
            img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                      local_DNAFISH_bg_correct, 0.01, 1)
        radial_distribution_relative_r_nuclear = \
            img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                      local_nuclear, 0.01, 1)

        radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
        radial_distribution_relative_r_DNAFISH_bg_correct_smooth = \
            dat.list_smooth(radial_distribution_relative_r_DNAFISH_bg_correct, 3)
        radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

        radial_subtract = \
            np.array(radial_distribution_relative_r_DNAFISH_bg_correct_smooth) - np.array(
                radial_distribution_relative_r_nuclear_smooth)
        radial_subtract_center = np.mean(radial_subtract[0:40])
        radial_subtract_edge = np.mean(radial_subtract[40:80])

        # angle distribution
        local_angle_map = img.angle_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
        angle_distribution_DNAFISH = img.radial_distribution_from_distance_map(local_nuclear_seg_convex,
                                                                               local_angle_map,
                                                                               local_DNAFISH, 1, 360)
        angle_distribution_DNAFISH_bg_correct = \
            img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map,
                                                      local_DNAFISH_bg_correct, 1, 360)
        angle_distribution_nuclear = img.radial_distribution_from_distance_map(local_nuclear_seg_convex,
                                                                               local_angle_map,
                                                                               local_nuclear, 1, 360)

        angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
        angle_distribution_DNAFISH_bg_correct_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH_bg_correct, 7)
        angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
        angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
            dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)
        angle_distribution_DNAFISH_bg_correct_smooth_centered = \
            dat.list_peak_center(angle_distribution_DNAFISH_bg_correct_smooth)

        data.loc[len(data.index)] = [i, data_z['FOV'][i], data_z['z'][i], data_z['z_min'][i], data_z['z_max'][i],
                                     data_z['z_ratio'][i], data_z['label_nuclear'][i],
                                     data_z['centroid_nuclear'][i], data_z['limit'][i], area_nuclear,
                                     mean_int_bg, mean_int_DNAFISH,
                                     mean_int_DNAFISH_norm,
                                     mean_int_nuclear, total_int_DNAFISH, total_int_DNAFISH_norm,
                                     total_int_nuclear, radial_distribution_relative_r_DNAFISH_smooth,
                                     radial_distribution_relative_r_DNAFISH_bg_correct_smooth,
                                     radial_distribution_relative_r_nuclear_smooth, radial_subtract_center,
                                     radial_subtract_edge, angle_distribution_DNAFISH_smooth_centered,
                                     angle_distribution_DNAFISH_bg_correct_smooth_centered,
                                     angle_distribution_nuclear_smooth_centered, total_area_ecDNA,
                                     total_area_ratio_ecDNA,
                                     total_int_ecDNA,
                                     total_int_ecDNA_norm,
                                     mean_int_ecDNA, mean_int_ecDNA_norm,
                                     area_ind_ecDNA, area_ratio_ind_ecDNA,
                                     mean_int_ind_ecDNA, mean_int_ind_ecDNA_norm, total_int_ind_ecDNA,
                                     total_int_ind_ecDNA_norm,
                                     cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                     cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                     cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                     cum_percentage_total_int_ind_ecDNA_norm, cum_percentage_total_int_norm_n_half,
                                     cum_area_ind_ecDNA, cum_area_n_half, cum_area_ratio_ind_ecDNA,
                                     cum_area_ratio_n_half, cum_int_ind_ecDNA, cum_int_n_half, cum_int_ind_ecDNA_norm,
                                     cum_int_norm_n_half]

    else:
        tif.imwrite('%s%s/%s/%s_DNAFISH_fov%s_z%s_i%s.tif' % (master_folder, sample_row, sample_num, sample, fov,
                                                              z_current, label_nuclear), local_DNAFISH)


data['max_area_ecDNA'] = [np.max(data['area_individual_ecDNA'][i]) for i in range(len(data))]
data['max_area_ratio_ecDNA'] = data['max_area_ecDNA']/data['area_nuclear']
data['n_ecDNA'] = [len(data['area_individual_ecDNA'][i]) for i in range(len(data))]
max_n_ecDNA = max(data['n_ecDNA'])
data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])
data['cum_int_ind_ecDNA_norm_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA_norm'])

data.to_csv('%s%s/%s/%s.txt' % (master_folder, sample_row, sample_num, sample), index=False, sep='\t')

print("DONE!")
