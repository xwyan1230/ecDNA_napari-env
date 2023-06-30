import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import os
import shared.dataframe as dat
import shared.objects as obj
import shared.math as mat
import cv2
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

row = 'C'
sample = 'C3'
batch = 1
total_fov = 4
dshape_factor = 0.0765
n_nuclear_convex_dilation = 4
local_size = 200
rmax = 100
start_fov = 0


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


df_align = pd.read_csv('%s/alignment/%s/%s_alignment_autocheck_%s.txt' % (master_folder, sample, sample, batch), na_values=['.'], sep='\t')

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (data_dir, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'FOV', 'label',
                             'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear',
                             'mean_int_red', 'mean_int_green',
                             'mean_int_DNAFISH', 'n_ecDNA',
                             'mean_int_ecDNA', 'area_ecDNA_individual', 'area_ratio_ecDNA_individual', 'mean_int_ecDNA_individual',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'dis_to_hub_area',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half',
                             'n_ecDNA_lst', 'total_area_DNAFISH_lst', 'total_area_ratio_DNAFISH_lst'])

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_4pos_s%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s%s' % (row, sample, fov)
    img_nuclear_bgc = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_%s_seg.tif" % (data_dir, sample, sample, batch, fov),
                                  plugin="tifffile")
    img_ecSeg = skio.imread("%s/seg/%s/seg_tif/%s_%s_%s_ecseg1.tif" % (data_dir, sample, sample, batch, fov),
                                plugin="tifffile")
    s = img_nuclear_bgc.shape[0]*dshape_factor
    topleft = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]
    img_before_GFP_cut = img_before_GFP.copy()[int(topleft[1]):int(topleft[1]+s)+5, int(topleft[0]):int(topleft[0]+s)+5]
    s1 = img_before_GFP_cut.shape[0]
    img_before_GFP_cut_resize = cv2.resize(img_before_GFP_cut, dsize=(int(s1 * 1/dshape_factor), int(s1 * 1/dshape_factor)),
                                           interpolation=cv2.INTER_AREA)[:img_nuclear_bgc.shape[0], :img_nuclear_bgc.shape[1]]
    img_before_mCherry_cut = img_before_mCherry.copy()[int(topleft[1]):int(topleft[1]+s)+5, int(topleft[0]):int(topleft[0]+s)+5]
    img_before_mCherry_cut_resize = cv2.resize(img_before_mCherry_cut,
                                           dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                           interpolation=cv2.INTER_AREA)[:img_nuclear_bgc.shape[0], :img_nuclear_bgc.shape[1]]

    """viewer = napari.Viewer()
    viewer.add_image(img_before_mCherry_cut_resize, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_before_GFP_cut_resize, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_seg, blending='additive', colormap='green', contrast_limits=[0, 20])
    napari.run()"""

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    # measure
    # get local images
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        label_nuclear = nuclear_props[i].label
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_red = img_before_mCherry_cut_resize.copy()[position[0]:position[1], position[2]:position[3]]
        local_green = img_before_GFP_cut_resize.copy()[position[0]:position[1], position[2]:position[3]]

        """viewer = napari.Viewer()
        viewer.add_image(local_red, blending='additive', colormap='red', contrast_limits=[0, 65535])
        viewer.add_image(local_green, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(local_nuclear_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
        napari.run()"""

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_red_props = regionprops(label(local_nuclear_seg), local_red)
        local_green_props = regionprops(label(local_nuclear_seg), local_green)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        perimeter_nuclear = local_nuclear_props[0].perimeter
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        mean_int_red = local_red_props[0].intensity_mean
        mean_int_green = local_green_props[0].intensity_mean
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

        n_ecDNA_lst = []
        total_area_DNAFISH_lst = []
        total_area_ratio_DNAFISH_lst = []
        for j in range(15):
            thresh = j * 5000
            local_DNAFISH_seg = np.zeros_like(local_DNAFISH)
            local_DNAFISH_seg[local_DNAFISH >= thresh] = 1
            local_DNAFISH_seg = obj.remove_small(local_DNAFISH_seg, 20)
            local_DNAFISH_seg[local_nuclear_seg == 0] = 0
            ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)
            n_ecDNA_lst.append(len(ecDNA_props))

            local_DNAFISH_seg = np.zeros_like(local_DNAFISH)
            local_DNAFISH_seg[local_DNAFISH >= thresh] = 1
            local_DNAFISH_seg[local_DNAFISH >= thresh + 5000] = 0
            local_DNAFISH_seg[local_nuclear_seg == 0] = 0
            total_area_DNAFISH_lst.append(local_DNAFISH_seg.sum())
            total_area_ratio_DNAFISH_lst.append(local_DNAFISH_seg.sum() / area_nuclear)

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
                                                      local_DNAFISH, 0.05, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                      local_nuclear, 0.05, 1)

        data.loc[len(data.index)] = [i, fov, label_nuclear,
                                     n_nuclear_convex_dilation,
                                     original_centroid_nuclear, area_nuclear, circ_nuclear, mean_int_nuclear,
                                     mean_int_red, mean_int_green,
                                     mean_int_DNAFISH, n_ecDNA,
                                     mean_int_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA, mean_int_ind_ecDNA,
                                     total_area_ecDNA, total_area_ratio_ecDNA,
                                     dis_to_hub_area_v2,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                     cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                     cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                     cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                     cum_area_ind_ecDNA, cum_area_n_half,
                                     cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half,
                                     cum_int_ind_ecDNA, cum_int_n_half,
                                     n_ecDNA_lst, total_area_DNAFISH_lst, total_area_ratio_DNAFISH_lst]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])

if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
data.to_csv('%s%s/%s_n%s_simplified_%s.txt' % (output_dir, sample, sample, n_nuclear_convex_dilation, batch), index=False, sep='\t')

print("DONE!")