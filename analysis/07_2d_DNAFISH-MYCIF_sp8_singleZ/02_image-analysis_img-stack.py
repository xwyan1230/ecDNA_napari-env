from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, binary_dilation, disk, binary_erosion
import pandas as pd
import numpy as np
import shared.image as img
import skimage.io as skio
import shared.image as ima
import shared.dataframe as dat
import shared.segmentation as seg
import shared.math as mat
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220916_double_funnel_test/"
sample = 'slide1_U2'
sample_folder = sample
save_path = master_folder
start_fov = 0
total_fov = 30
version = 1

local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV',
                             'centroid_nuclear', 'area_nuclear',
                             'bg_int_nuclear', 'mean_int_nuclear', 'total_int_nuclear',
                             'bg_int_IF', 'mean_int_IF', 'total_int_IF',
                             'bg_int_DNAFISH', 'mean_int_DNAFISH', 'total_int_DNAFISH',
                             'n_ecDNA',
                             'mean_int_ecDNA', 'total_int_ecDNA', 'mean_int_ind_ecDNA', 'total_int_ind_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                             'max_area_ecDNA', 'max_area_ratio_ecDNA',
                             'g', 'g_value',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_center', 'radial_edge',
                             'relative_r_area', 'relative_r_int',
                             'angle_curve_nuclear', 'angle_curve_DNAFISH', 'angle_value',
                             'dis_to_hub_area', 'dis_to_hub_int',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half'])

# load images
file_prefix = '20220916_double_funnel_test_%s_RAW' % sample
img_nuclear_stack = skio.imread("%s%s/%s_ch00.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch01.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_IF_stack = skio.imread("%s%s/%s_ch02.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_nuclear_seg_stack = skio.imread("%s%s/%s_seg.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")
img_DNAFISH_seg_stack = skio.imread("%s%s/%s_ecSeg.tif" % (master_folder, sample_folder, file_prefix), plugin="tifffile")

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, fov %s/%s" % (sample, fov+1, total_fov))

    # load images
    img_nuclear = img_nuclear_stack[fov]
    img_DNAFISH = img_DNAFISH_stack[fov]
    img_IF = img_IF_stack[fov]
    img_nuclear_seg = img_nuclear_seg_stack[fov]
    img_DNAFISH_seg = img_DNAFISH_seg_stack[fov]
    """for r in range(5):
        img_DNAFISH_seg = binary_dilation(img_DNAFISH_seg)
    for r in range(5):
        img_DNAFISH_seg = binary_erosion(img_DNAFISH_seg)"""

    # background measurement and background correction
    _, sample_img_nuclear = seg.get_bg_img(img_nuclear)
    _, sample_img_DNAFISH = seg.get_bg_img(img_DNAFISH)
    _, sample_img_IF = seg.get_bg_img(img_IF)

    sample_seg = np.zeros_like(img_nuclear)
    sample_seg[sample_img_nuclear == 1] = 1
    sample_seg[sample_img_DNAFISH == 1] = 1
    sample_seg[sample_img_IF == 1] = 1
    sample_seg = binary_dilation(sample_seg, disk(50))
    bg = np.ones_like(sample_seg)
    bg[sample_seg == 1] = 0

    bg_int_nuclear = np.sum(bg * img_nuclear) / np.sum(bg)
    bg_int_DNAFISH = np.sum(bg * img_DNAFISH) / np.sum(bg)
    bg_int_IF = np.sum(bg * img_IF) / np.sum(bg)
    img_nuclear_bgc = img_nuclear
    img_DNAFISH_bgc = img_DNAFISH
    img_IF_bgc = img_IF

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
        local_IF = img_IF_bgc.copy()
        local_IF = local_IF[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_DNAFISH_seg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        local_IF_props = regionprops(label(local_nuclear_seg), local_IF)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        total_int_nuclear = area_nuclear * mean_int_nuclear
        mean_int_IF = local_IF_props[0].intensity_mean
        total_int_IF = area_nuclear * mean_int_IF
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        total_int_DNAFISH = area_nuclear * mean_int_DNAFISH

        # ecDNA measurements
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

        # distance from the hub
        dis_to_hub_area = 0
        dis_to_hub_int = 0

        if n_ecDNA == 0:
            dis_to_hub_area = -1
            dis_to_hub_int = -1
        elif n_ecDNA > 1:
            ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'total_int': total_int_ind_ecDNA,
                                      'centroid': centroid_ind_ecDNA})

            ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
                                                               ignore_index=True)
            ind_ecDNA_sort_area['dis'] = \
                [((ind_ecDNA_sort_area['centroid'][i][0] - ind_ecDNA_sort_area['centroid'][0][0]) ** 2 +
                  (ind_ecDNA_sort_area['centroid'][i][1] - ind_ecDNA_sort_area['centroid'][0][1]) ** 2) ** 0.5
                 for i in range(len(ind_ecDNA_sort_area))]

            for ind in range(n_ecDNA - 1):
                dis_to_hub_area = dis_to_hub_area + ind_ecDNA_sort_area['area'][ind + 1] / total_area_ecDNA * \
                                  ind_ecDNA_sort_area['dis'][ind + 1]

            ind_ecDNA_sort_int = ind_ecDNA.copy().sort_values(by='total_int', axis=0, ascending=False, inplace=False,
                                                              ignore_index=True)
            ind_ecDNA_sort_int['dis'] = \
                [((ind_ecDNA_sort_int['centroid'][i][0] - ind_ecDNA_sort_int['centroid'][0][0]) ** 2 +
                  (ind_ecDNA_sort_int['centroid'][i][1] - ind_ecDNA_sort_int['centroid'][0][1]) ** 2) ** 0.5
                 for i in range(len(ind_ecDNA_sort_int))]

            for ind in range(n_ecDNA - 1):
                dis_to_hub_int = dis_to_hub_int + ind_ecDNA_sort_int['total_int'][ind + 1] / total_int_ecDNA * \
                                 ind_ecDNA_sort_int['dis'][ind + 1]

        # auto-correlation
        _, r, g, dg = mat.auto_correlation(local_DNAFISH, local_nuclear_seg, rmax)
        g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = img.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            img.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH, 0.01, 1)
        radial_distribution_relative_r_nuclear = \
            img.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_nuclear, 0.01, 1)

        radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
        radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

        radial_subtract = np.array(radial_distribution_relative_r_DNAFISH_smooth) - \
                          np.array(radial_distribution_relative_r_nuclear_smooth)
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
        local_angle_map = img.angle_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        angle_distribution_DNAFISH = \
            img.radial_distribution_from_distance_map(local_nuclear_seg, local_angle_map, local_DNAFISH, 1, 360)
        angle_distribution_nuclear = \
            img.radial_distribution_from_distance_map(local_nuclear_seg, local_angle_map, local_nuclear, 1, 360)

        angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
        angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
        angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
            dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)
        angle_value = angle_distribution_DNAFISH_smooth_centered[179]

        data.loc[len(data.index)] = [i, fov,
                                     original_centroid_nuclear, area_nuclear,
                                     bg_int_nuclear, mean_int_nuclear, total_int_nuclear,
                                     bg_int_IF, mean_int_IF, total_int_IF,
                                     bg_int_DNAFISH, mean_int_DNAFISH, total_int_DNAFISH,
                                     n_ecDNA,
                                     mean_int_ecDNA, total_int_ecDNA, mean_int_ind_ecDNA, total_int_ind_ecDNA,
                                     total_area_ecDNA, total_area_ratio_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA,
                                     max_area_ecDNA, max_area_ratio_ecDNA,
                                     g, g_value,
                                     radial_distribution_relative_r_nuclear_smooth,
                                     radial_distribution_relative_r_DNAFISH_smooth,
                                     radial_subtract_center, radial_subtract_edge,
                                     relative_r_area, relative_r_int,
                                     angle_distribution_nuclear_smooth_centered,
                                     angle_distribution_DNAFISH_smooth_centered, angle_value,
                                     dis_to_hub_area, dis_to_hub_int,
                                     cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                     cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                     cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                     cum_area_ind_ecDNA, cum_area_n_half,
                                     cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half,
                                     cum_int_ind_ecDNA, cum_int_n_half]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])

data.to_csv('%s%s_v%s.txt' % (master_folder, sample, version), index=False, sep='\t')

print("DONE!")



