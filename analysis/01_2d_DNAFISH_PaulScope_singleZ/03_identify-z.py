import pandas as pd
import skimage.io as skio
import shared.dataframe as dat
import shared.image as img
import numpy as np
from skimage.filters import threshold_li
import napari

# INPUT PARAMETERS
# file info
# master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220816_Natasha_ColoDM_reimage/"
master_folder = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220825_POLR3D/20220825_POLR3Dtest/'
sample = 'Control_2'
master_path = '%s%s/' % (master_folder, sample)
total_fov = 10
start_fov = 1
# centroid searching
nuclear_centroid_searching_range = 25  # pixel
local_size = 100

# LOAD CENTROIDS FILE
data_c = pd.read_csv('%s%s_centroids.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data_c['centroid_nuclear'] = [dat.str_to_float(data_c['centroid_nuclear'][i]) for i in range(len(data_c))]
data_c['centroid_x'] = [data_c['centroid_nuclear'][i][0] for i in range(len(data_c))]
data_c['centroid_y'] = [data_c['centroid_nuclear'][i][1] for i in range(len(data_c))]

data = pd.DataFrame(columns=['nuclear', 'FOV', 'z_min', 'z_max', 'z', 'label_nuclear', 'centroid_nuclear', 'limit',
                             'z_ratio'])

n_nuclear = pd.DataFrame(columns=['sample', 'FOV', 'total', 'fov_skipped', 'surrounding_z', 'duplicate_z',
                                  'discontinue_z', 'fewer_z', 'filtered_total'])

# IMAGING ANALYSIS
for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, start nuclear analysis FOV %s/%s" % (sample, f+1, total_fov))
    file_prefix = "%s_Position %s_RAW" % (sample, fov)
    duplicate_z = 0
    discontinue_z = 0
    fewer_z = 0
    filtered_total = 0
    surrounding_z = 0
    # load images
    """im_z_stack_nuclear = img.img_to_int(skio.imread("%s%s/%s_ch00.tif" % (master_path, raw_folder, file_prefix),
                                                        plugin="tifffile"))

        im_z_stack_DNAFISH = img.img_to_int(skio.imread("%s%s/%s_ch01.tif" % (master_path, raw_folder, file_prefix),
                                                        plugin="tifffile"))"""
    im_z_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s_ch01.tif" % (master_path, file_prefix), plugin="tifffile")

    total_z = im_z_stack_nuclear.shape[0]
    im_z_stack_seg_convex = skio.imread("%s%s_seg.tif" % (master_path, file_prefix), plugin="tifffile")
    # identify z for given fov
    data_c_fov = data_c[data_c['FOV'] == fov]
    z_analyze = int(total_z/2)
    data_c_fov_z = data_c_fov[data_c_fov['z'] == z_analyze]
    n_nuclear_fov = len(data_c_fov[data_c_fov['z'] == z_analyze])

    # checking for drying situation
    threshold_li_cutoff = threshold_li(im_z_stack_nuclear[z_analyze])
    bg_li_seg = im_z_stack_nuclear[z_analyze] < threshold_li_cutoff
    fov_bg = im_z_stack_nuclear[z_analyze] * bg_li_seg
    mean_fov_bg = np.sum(fov_bg) / np.sum(bg_li_seg)
    if mean_fov_bg > 0:
        # identify each nucleus
        for i in range(n_nuclear_fov):
            print("Analyzing nucleus %s/%s" % (i + 1, n_nuclear_fov))
            data_nucleus = data_c_fov[(data_c_fov['centroid_x'] > data_c_fov_z['centroid_x'].tolist()[i] -
                                       nuclear_centroid_searching_range)
                                      & (data_c_fov['centroid_x'] < data_c_fov_z['centroid_x'].tolist()[i] +
                                         nuclear_centroid_searching_range)
                                      & (data_c_fov['centroid_y'] > data_c_fov_z['centroid_y'].tolist()[i] -
                                         nuclear_centroid_searching_range)
                                      & (data_c_fov['centroid_y'] < data_c_fov_z['centroid_y'].tolist()[i] +
                                         nuclear_centroid_searching_range)].copy()
            z_min = data_nucleus['z'].tolist()[0]
            z_max = data_nucleus['z'].tolist()[-1]
            if (len(data_nucleus) == len(set(data_nucleus['z'].tolist()))) & (len(data_nucleus) >= 3):
                if (len(data_nucleus) == z_max - z_min + 1) | (len(data_nucleus) >= 10):
                    filtered_total = filtered_total + 1
                    limit = int(data_nucleus['mean_intensity_MYC_DNAFISH_in_nucleus'].tolist()[0] * 2)

                    mean_int_DNAFISH_over_limit = []
                    for z in range(len(data_nucleus)):
                        # find nucleus for given z
                        z_current = data_nucleus['z'].tolist()[z]
                        label_nuclear = data_nucleus['label_nuclear'].tolist()[z]
                        original_centroid_nuclear = data_nucleus['centroid_nuclear'].tolist()[z]
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
                        local_DNAFISH[local_nuclear_seg_convex == 0] = 0
                        # calculate mean_int of DNAFISH signal higher than limit
                        local_DNAFISH[local_DNAFISH <= limit] = 0
                        DNAFISH_over_limit = local_DNAFISH.ravel()
                        DNAFISH_over_limit = [k for k in DNAFISH_over_limit if k != 0]
                        if len(DNAFISH_over_limit) == 0:
                            mean_int_DNAFISH_over_limit.append(0)
                        else:
                            mean_int_DNAFISH_over_limit.append(np.mean(DNAFISH_over_limit))

                    data_nucleus.loc[:, ['mean_intensity_over_limit']] = mean_int_DNAFISH_over_limit
                    data_nucleus_sort = data_nucleus.copy().sort_values(by='mean_intensity_over_limit', ascending=False)
                    data_nucleus_sort.reset_index(drop=True, inplace=True)
                    z_signal = data_nucleus_sort['z'].tolist()[0]
                    z_ratio = (z_signal - z_min) / (z_max - z_min)
                    if (z_signal - 1 in data_nucleus['z'].tolist()) & (z_signal - 2 in data_nucleus['z'].tolist()) & \
                            (z_signal + 1 in data_nucleus['z'].tolist()) & (z_signal + 2 in data_nucleus['z'].tolist()):
                        if z_signal != z_min:
                            data.loc[len(data.index)] = [i, fov, z_min, z_max, z_signal,
                                                         data_nucleus_sort['label_nuclear'].tolist()[0],
                                                         data_nucleus_sort['centroid_nuclear'].tolist()[0], limit,
                                                         z_ratio]
                        else:
                            mid = int(len(data_nucleus) / 2)
                            data.loc[len(data.index)] = [i, fov, z_min, z_max, data_nucleus['z'].tolist()[mid],
                                                         data_nucleus['label_nuclear'].tolist()[mid],
                                                         data_nucleus['centroid_nuclear'].tolist()[mid], limit, z_ratio]
                    else:
                        surrounding_z = surrounding_z + 1
                        filtered_total = filtered_total - 1
                else:
                    discontinue_z = discontinue_z + 1
            elif len(data_nucleus) != len(set(data_nucleus['z'].tolist())):
                duplicate_z = duplicate_z + 1
            elif len(data_nucleus) < 3:
                fewer_z = fewer_z + 1
        n_nuclear.loc[len(n_nuclear)] = [sample, fov, n_nuclear_fov, 0, surrounding_z, duplicate_z, discontinue_z,
                                         fewer_z, filtered_total]

    else:
        n_nuclear.loc[len(n_nuclear)] = [sample, fov, n_nuclear_fov, 1, surrounding_z, duplicate_z, discontinue_z,
                                         fewer_z, filtered_total]

data.to_csv('%s%s_z.txt' % (master_folder, sample), index=False, sep='\t')
n_nuclear.to_csv('%s%s_n_nuclear.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")
