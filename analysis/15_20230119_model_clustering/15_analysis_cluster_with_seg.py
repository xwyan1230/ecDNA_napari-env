import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.math as mat
import shared.dataframe as dat
from skimage.morphology import medial_axis
import os
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%ssimulated_data/dataset2_different-r-nuclear_random-radial_different-ac_cp-100/" % master_folder
output_dir = "%stxt/dataset2/" % master_folder

sample = '5000_5'

local_size = 150
rmax = 75

sample_folder1 = '%s%s/FISH_seg_tif/' % (data_dir, sample)
sample_folder2 = '%s%s/nuclei_seg_tif/' % (data_dir, sample)
FISH_imgs = [x for x in os.listdir(sample_folder1)]
if '.DS_Store' in FISH_imgs:
    FISH_imgs.remove('.DS_Store')

data = pd.DataFrame(columns=['FOV', 'coefficient', 'crange', 'copy_num',
                             'area_nuclear',
                             'n_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'dis_to_hub_area_v2', 'g',
                             'percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA', 'radial_curve_DNAFISH'])

for i in range(len(FISH_imgs)):
    im_FISH = skio.imread("%s%s" % (sample_folder1, FISH_imgs[i]), plugin="tifffile")
    im_seg = skio.imread("%sim_seg_%s" % (sample_folder2, FISH_imgs[i][8:]), plugin="tifffile")
    fov = FISH_imgs[i].split('_')[2]
    coefficient = sample.split('_')[0]
    crange = sample.split('_')[1]
    copy_num = FISH_imgs[i].split('_')[3].split('.')[0][2:]
    print('sample: %s, cell: %s, copy number: %s' % (sample, i, copy_num))

    # measure
    area_nuclear = np.sum(im_seg)
    ecDNA_props = regionprops(label(im_FISH))
    n_ecDNA = len(ecDNA_props)
    centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
    area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
    area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
    total_area_ecDNA = sum(area_ind_ecDNA)
    total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

    # distance from hub v2
    dis_to_hub_area_v2 = 0

    if n_ecDNA == 0:
        dis_to_hub_area_v2 = 0
    elif n_ecDNA > 1:
        ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'centroid': centroid_ind_ecDNA})
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
                dis_temp = dis_temp + current_ecDNA['area'][ind + 1] * current_ecDNA['dis'][
                    ind + 1] / total_area_ecDNA
            dis.append(dis_temp)
        ind_ecDNA_sort_area['dis'] = dis

        for ind in range(n_ecDNA):
            dis_to_hub_area_v2 = dis_to_hub_area_v2 + ind_ecDNA_sort_area['area'][ind] * \
                                 ind_ecDNA_sort_area['dis'][ind] / total_area_ecDNA

    # auto-correlation
    _, r, g, dg = mat.auto_correlation(im_FISH, im_seg, rmax)

    # percentage/cumulative
    percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True)) / total_area_ecDNA)
    cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
    cum_area_ind_ecDNA = dat.list_sum(area_ind_ecDNA)

    # radial distribution
    local_nuclear_centroid = (local_size, local_size)
    _, local_edge_distance_map = medial_axis(im_seg, return_distance=True)
    local_centroid_distance_map = ima.distance_map_from_point(im_seg, local_nuclear_centroid)
    local_centroid_distance_map[im_seg == 0] = 0
    local_edge_distance_map[im_seg == 0] = -1
    local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)
    radial_distribution_relative_r_DNAFISH = \
        ima.radial_distribution_from_distance_map(im_seg, local_relative_r_map, im_FISH, 0.05, 1)

    data.loc[len(data.index)] = [fov, coefficient, crange, copy_num,
                                 area_nuclear,
                                 n_ecDNA, total_area_ecDNA, total_area_ratio_ecDNA,
                                 dis_to_hub_area_v2, g,
                                 cum_percentage_area_ind_ecDNA, cum_area_ind_ecDNA, radial_distribution_relative_r_DNAFISH]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data.to_csv('%s%s_cluster.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")
