import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis, erosion
import shared.objects as obj
import shared.image as ima
import shared.dataframe as dat
from skimage import segmentation
import numpy as np
import matplotlib
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import skimage.io as skio
import shared.math as mat
import math
import shared.display as dis
import pandas as pd
import tifffile as tif
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/E/"
sample_lst = ['HZ', 'FM', 'KW', 'FD', 'ME', 'CK', 'AO', 'IK', 'CR', 'BD', 'HD', 'NP', 'BQ', 'CO', 'LM', 'GP',
              'MM', 'LY', 'JZ', 'NB', 'LT', 'MK', 'DA', 'LX', 'GR', 'NQ', 'JF', 'DT']
group = 'E'
save_path = master_folder
save_folder = master_folder
version = 1

# default setting
hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]
line_color_4 = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70), (0.95, 0.85, 0)]

# other parameters
local_size = 100
rmax = 80

data = pd.DataFrame(columns=['sample', 'nuclear',
                             'centroid_nuclear', 'area_nuclear', 'normalized_r',
                             'bg_int_nuclear', 'mean_int_nuclear', 'total_int_nuclear',
                             'bg_int_DNAFISH', 'mean_int_DNAFISH', 'total_int_DNAFISH',
                             'n_ecDNA',
                             'mean_int_ecDNA', 'total_int_ecDNA', 'mean_int_ind_ecDNA', 'total_int_ind_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                             'max_area_ecDNA', 'max_area_ratio_ecDNA',
                             'g', 'g_value',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized',
                             'radial_curve_mask',
                             'relative_r_area', 'relative_r_int',
                             'angle_curve_nuclear', 'angle_curve_DNAFISH', 'angle_value',
                             'dis_to_hub_area', 'dis_to_hub_int', 'dis_to_hub_area_v2', 'normalized_cluster',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half'])

for sample in sample_lst:
    # load images
    img_nuclear = skio.imread("%s%s_b.BMP" % (master_folder, sample))[:, :, 2]
    img_DNAFISH = skio.imread("%s%s_g.BMP" % (master_folder, sample))[:, :, 1]
    img_nuclear_seg = skio.imread("%s%s_seg.tif" % (master_folder, sample), plugin="tifffile")
    img_DNAFISH_seg = skio.imread("%s%s_ecSeg.tif" % (master_folder, sample), plugin="tifffile")

    # bg_correction
    _, sample_img_nuclear = seg.get_bg_img(img_nuclear)
    _, sample_img_DNAFISH = seg.get_bg_img(img_DNAFISH)

    sample_seg = np.zeros_like(img_nuclear)
    sample_seg[sample_img_nuclear == 1] = 1
    sample_seg[sample_img_DNAFISH == 1] = 1
    sample_seg = binary_dilation(sample_seg, disk(50))
    bg = np.ones_like(sample_seg)
    bg[sample_seg == 1] = 0

    bg_int_nuclear = np.sum(bg * img_nuclear) / np.sum(bg)
    bg_int_DNAFISH = np.sum(bg * img_DNAFISH) / np.sum(bg)
    img_nuclear_bgc = img_nuclear
    img_DNAFISH_bgc = img_DNAFISH

    img_nuclear_bgc = img_nuclear.astype(float) - np.ones_like(img_nuclear) * bg_int_nuclear
    img_nuclear_bgc[img_nuclear_bgc < 0] = 0
    img_DNAFISH_bgc = img_DNAFISH.astype(float) - np.ones_like(img_DNAFISH) * bg_int_DNAFISH
    img_DNAFISH_bgc[img_DNAFISH_bgc < 0] = 0

    # get local images
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, nuclear %s/%s" % (sample, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_DNAFISH_seg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        total_int_nuclear = area_nuclear * mean_int_nuclear
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        total_int_DNAFISH = area_nuclear * mean_int_DNAFISH
        normalized_r = (area_nuclear/math.pi)**0.5

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

        # distance from hub v2
        dis_to_hub_area_v2 = 0

        if n_ecDNA == 0:
            dis_to_hub_area_v2 = 0
        elif n_ecDNA > 1:
            ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'total_int': total_int_ind_ecDNA,
                                      'centroid': centroid_ind_ecDNA})

            ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
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
                    dis_temp = dis_temp + current_ecDNA['area'][ind + 1] * current_ecDNA['dis'][ind+1] / total_area_ecDNA
                dis.append(dis_temp)
            ind_ecDNA_sort_area['dis'] = dis

            for ind in range(n_ecDNA):
                dis_to_hub_area_v2 = dis_to_hub_area_v2 + ind_ecDNA_sort_area['area'][ind] * ind_ecDNA_sort_area['dis'][ind] / total_area_ecDNA

        normalized_cluster = (dis_to_hub_area_v2/normalized_r)/(total_area_ratio_ecDNA**0.5)

        # auto-correlation
        _, r, g, dg = mat.auto_correlation(local_DNAFISH, local_nuclear_seg, rmax)
        g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_DNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH, 0.1, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_nuclear, 0.1, 1)

        # radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
        # radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

        radial_distribution_relative_r_normalized = \
            list(np.array(radial_distribution_relative_r_DNAFISH)/np.array(radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_normalized = dat.nan_replace(radial_distribution_relative_r_normalized)
        if ~np.isnan(radial_distribution_relative_r_normalized).any():
            # radial distribution of mask
            radial_distribution_relative_r_mask = \
                ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH_seg, 0.1, 1)

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
            local_angle_map = ima.angle_map_from_point(local_nuclear_seg, local_nuclear_centroid)
            angle_distribution_DNAFISH = \
                ima.radial_distribution_from_distance_map(local_nuclear_seg, local_angle_map, local_DNAFISH, 1, 360)
            angle_distribution_nuclear = \
                ima.radial_distribution_from_distance_map(local_nuclear_seg, local_angle_map, local_nuclear, 1, 360)

            angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
            angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
            angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
                dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)
            angle_value = angle_distribution_DNAFISH_smooth_centered[179]

            data.loc[len(data.index)] = [sample, i,
                                         original_centroid_nuclear, area_nuclear, normalized_r,
                                         bg_int_nuclear, mean_int_nuclear, total_int_nuclear,
                                         bg_int_DNAFISH, mean_int_DNAFISH, total_int_DNAFISH,
                                         n_ecDNA,
                                         mean_int_ecDNA, total_int_ecDNA, mean_int_ind_ecDNA, total_int_ind_ecDNA,
                                         total_area_ecDNA, total_area_ratio_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA,
                                         max_area_ecDNA, max_area_ratio_ecDNA,
                                         g, g_value,
                                         radial_distribution_relative_r_nuclear,
                                         radial_distribution_relative_r_DNAFISH,
                                         radial_distribution_relative_r_normalized,
                                         radial_distribution_relative_r_mask,
                                         relative_r_area, relative_r_int,
                                         angle_distribution_nuclear_smooth_centered,
                                         angle_distribution_DNAFISH_smooth_centered, angle_value,
                                         dis_to_hub_area, dis_to_hub_int, dis_to_hub_area_v2, normalized_cluster,
                                         cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                         cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                         cum_percentage_total_int_ind_ecDNA, cum_percentage_total_int_n_half,
                                         cum_area_ind_ecDNA, cum_area_n_half,
                                         cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half,
                                         cum_int_ind_ecDNA, cum_int_n_half]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])
data['cum_int_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_int_ind_ecDNA'])
data_drop = data[data['area_nuclear'] > 100].copy().reset_index(drop=True)

data_drop.to_csv('%s%s_v%s.txt' % (master_folder, group, version), index=False, sep='\t')

data = data_drop.copy()
# heatmap
data_heatmap = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
sample_updated = []
sample_n = []
for sample in sample_lst:
    data_sample = data[data['sample'] == sample].copy().reset_index(drop=True)
    if len(data_sample) >= 10:
        sample_updated.append(sample)
        sample_n.append(len(data_sample))
        data_radial = pd.DataFrame()
        data_radial['0.05'] = [data_sample['radial_curve_normalized'][i][0] for i in range(len(data_sample))]
        data_radial['0.15'] = [data_sample['radial_curve_normalized'][i][1] for i in range(len(data_sample))]
        data_radial['0.25'] = [data_sample['radial_curve_normalized'][i][2] for i in range(len(data_sample))]
        data_radial['0.35'] = [data_sample['radial_curve_normalized'][i][3] for i in range(len(data_sample))]
        data_radial['0.45'] = [data_sample['radial_curve_normalized'][i][4] for i in range(len(data_sample))]
        data_radial['0.55'] = [data_sample['radial_curve_normalized'][i][5] for i in range(len(data_sample))]
        data_radial['0.65'] = [data_sample['radial_curve_normalized'][i][6] for i in range(len(data_sample))]
        data_radial['0.75'] = [data_sample['radial_curve_normalized'][i][7] for i in range(len(data_sample))]
        data_radial['0.85'] = [data_sample['radial_curve_normalized'][i][8] for i in range(len(data_sample))]
        data_radial['0.95'] = [data_sample['radial_curve_normalized'][i][9] for i in range(len(data_sample))]

        percentage_higher = [len(data_radial[data_radial['0.05'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.15'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.25'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.35'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.45'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.55'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.65'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.75'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.85'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.95'] > 1])/len(data_radial)]
        data_heatmap.loc[len(data_heatmap.index)] = percentage_higher

print(sample_updated)
print(sample_n)

# hierarchical Clustering
# AHC (agglomerative, bottom-up)
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.

plt.figure(figsize=(12, 9))
clusters = shc.linkage(data_heatmap, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/%s_ahc.pdf' % (save_folder, group))
plt.close()

nodes = R['ivl']
sample_nodes = [sample_updated[int(i)] for i in nodes]

data_heatmap_sort = pd.DataFrame(columns=data_heatmap.columns)
for i in range(len(data_heatmap)):
    data_heatmap_sort.loc[len(data_heatmap_sort.index)] = \
        data_heatmap.iloc[int(dat.list_invert(nodes)[i])]

# heat map
plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap.pdf' % (save_folder, group))
plt.close()

plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap_sort, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='coolwarm',
                  yticklabels=False)
plt.savefig('%s/%s_heatmap_sort.pdf' % (save_folder, group))
plt.close()

print("DONE!")