import pandas as pd
import skimage.io as skio
import shared.dataframe as dat
import shared.image as img
import numpy as np
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = '72hr_100nMTHZ'
raw_folder = '01_raw'
seg_folder = '02_seg'
total_fov = 16
total_z = 18
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
# centroid searching
nuclear_centroid_searching_range = 25  # pixel
local_size = 100

# LOAD CENTROIDS FILE
data_c = pd.read_csv('%s%s/%s_centroids.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_c['centroid_nuclear'] = [dat.str_to_float(data_c['centroid_nuclear'][i]) for i in range(len(data_c))]
data_c['centroid_x'] = [data_c['centroid_nuclear'][i][0] for i in range(len(data_c))]
data_c['centroid_y'] = [data_c['centroid_nuclear'][i][1] for i in range(len(data_c))]

data = pd.DataFrame(columns=['nuclear',
                             'FOV',
                             'z_min',
                             'z_max',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'limit'])

# IMAGING ANALYSIS
for fov in range(total_fov):
    print("Start nuclear analysis FOV %s/%s" % (fov + 1, total_fov))
    # load images
    im_z_stack_nuclear = skio.imread("%s%s/%s/%s_RAW_ch00_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s_RAW_ch01_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    im_z_stack_seg_convex = skio.imread("%s%s/%s/%s_seg_fov%s.tif" % (master_folder, sample, seg_folder, sample, fov),
                                        plugin="tifffile")
    # identify z for given fov
    data_c_fov = data_c[data_c['FOV'] == fov]
    z_analyze = int(total_z/2)
    data_c_fov_z = data_c_fov[data_c_fov['z'] == z_analyze]
    n_nuclear_fov = len(data_c_fov[data_c_fov['z'] == z_analyze])

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
        if (len(data_nucleus) == len(set(data_nucleus['z'].tolist()))) & (len(data_nucleus) >= 3):
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
            z_min = data_nucleus['z'].tolist()[0]
            z_max = data_nucleus['z'].tolist()[-1]
            data_nucleus = data_nucleus.sort_values(by='mean_intensity_over_limit', ascending=False)
            data_nucleus.reset_index(drop=True, inplace=True)

            data.loc[len(data.index)] = [i, fov, z_min, z_max, data_nucleus['z'].tolist()[0],
                                         data_nucleus['label_nuclear'].tolist()[0],
                                         data_nucleus['centroid_nuclear'].tolist()[0], limit]

data['z_ratio'] = (data['z'] - data['z_min'])/(data['z_max']-data['z_min'])

data.to_csv('%s%s/%s_z.txt' % (master_folder, sample, sample), index=False, sep='\t')

print("DONE!")
