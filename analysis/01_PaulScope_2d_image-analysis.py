import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import numpy as np
import shared.image as img
import napari
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
n_nuclear_convex_dilation = 3
nuclear_centroid_searching_range = 25  # pixel
local_size = 100

# SET UP PARAMETERS
total_fov = 1
total_z = 22

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s/%s_RAW_ch00_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s/%s_RAW_ch01_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_seg = skio.imread("%s%s/%s_seg_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")

data_z = pd.read_csv('%s%s/%s_fov3_z.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

data = pd.DataFrame(columns=['nuclear',
                             'FOV',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'limit',
                             'area_nuclear',
                             'total_int_DNAFISH',
                             'total_int_nuclear',
                             'radial_curve_DNAFISH',
                             'radial_curve_nuclear',
                             'radial_center',
                             'radial_edge',
                             'angle_curve_DNAFISH',
                             'angle_curve_nuclear',
                             'total_area_ecDNA',
                             'total_int_ecDNA',
                             'percentage_area_curve_ecDNA',
                             'percentage_int_curve_ecDNA'
                             ])

for i in range(len(data_z)):
    print("Analyzing nucleus %s/%s" % (i+1, len(data_z)))
    z = data_z['z'][i]
    label_nuclear = data_z['label_nuclear'][i]
    original_centroid_nuclear = data_z['centroid_nuclear'][i]

    img_nuclear_seg_convex = im_stack_seg[z]
    img_nuclear = im_stack_nuclear[z]
    img_DNAFISH = im_stack_DNAFISH[z]

    position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
    local_nuclear = img_nuclear.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]

    local_nuclear_props = regionprops(label(local_nuclear_seg_convex), local_nuclear)
    area_nuclear = local_nuclear_props[0].area
    total_int_nuclear = area_nuclear * local_nuclear_props[0].intensity_mean
    local_DNAFISH_props = regionprops(label(local_nuclear_seg_convex), local_DNAFISH)
    total_int_DNAFISH = area_nuclear * local_DNAFISH_props[0].intensity_mean

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
    radial_distribution_relative_r_nuclear = \
        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                  local_nuclear, 0.01, 1)

    radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
    radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

    radial_subtract = \
        np.array(radial_distribution_relative_r_DNAFISH_smooth)-np.array(radial_distribution_relative_r_nuclear_smooth)
    radial_subtract_center = np.mean(radial_subtract[0:40])
    radial_subtract_edge = np.mean(radial_subtract[40:80])

    # angle distribution
    local_angle_map = img.angle_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
    angle_distribution_DNAFISH = img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map,
                                                                           local_DNAFISH, 1, 360)
    angle_distribution_nuclear = img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map,
                                                                           local_nuclear, 1, 360)

    angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
    angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
    angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
        dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)

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
    DNAFISH_seg, _ = seg.find_organelle(local_DNAFISH, 'na', extreme_val=500, bg_val=data_z['limit'][i]*0.8, min_size=0,
                                        max_size=10000)
    DNAFISH_seg[local_nuclear_seg_convex == 0] = 0
    ecDNA_seg = img_dot_seg.copy()
    ecDNA_seg[DNAFISH_seg == 1] = 1

    ecDNA_props = regionprops(label(ecDNA_seg), local_DNAFISH)
    n_ecDNA = len(ecDNA_props)
    area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
    total_area_ecDNA = sum(area_ind_ecDNA)
    mean_int_ind_ecDNA = [ecDNA_props[i].intensity_mean for i in range(len(ecDNA_props))]
    total_int_ind_ecDNA = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA))
    total_int_ecDNA = sum(total_int_ind_ecDNA)

    percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True))/total_area_ecDNA)
    percentage_total_int_ind_ecDNA = list(np.array(sorted(total_int_ind_ecDNA, reverse=True))/total_int_ecDNA)
    cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
    cum_percentage_total_int_ind_ecDNA = dat.list_sum(percentage_total_int_ind_ecDNA)

    """viewer = napari.view_image(local_DNAFISH, blending='additive', colormap='green')
    viewer.add_image(ecDNA_seg, blending='additive', colormap='red')
    napari.run()
    ax = plt.plot(cum_percentage_area_ind_ecDNA, c='g')
    plt.plot(cum_percentage_total_int_ind_ecDNA, c='b')
    plt.show()"""

    data.loc[len(data.index)] = [i, data_z['FOV'][i], data_z['z'][i], data_z['label_nuclear'][i],
                                 data_z['centroid_nuclear'][i], data_z['limit'][i], area_nuclear, total_int_DNAFISH,
                                 total_int_nuclear, radial_distribution_relative_r_DNAFISH_smooth,
                                 radial_distribution_relative_r_nuclear_smooth, radial_subtract_center,
                                 radial_subtract_edge, angle_distribution_DNAFISH_smooth_centered,
                                 angle_distribution_nuclear_smooth_centered, total_area_ecDNA, total_int_ecDNA,
                                 cum_percentage_area_ind_ecDNA, cum_percentage_total_int_ind_ecDNA]

data.to_csv('%s%s/%s_fov3.txt' % (master_folder, sample, sample), index=False, sep='\t')
print("DONE!")
