import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


sample = 'PC9'
# fov_lst = [fov_to_str(i) for i in np.arange(2, 37, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 7, 1)]  # Colo320DM
# fov_lst = [fov_to_str(i) for i in np.arange(1, 31, 1)]  # Colo320HSR
# fov_lst = ['MycI2_'+fov_to_str(i) for i in np.arange(1, 17, 1)]  # HCT116
# fov_lst = list(np.arange(1, 9, 1)) + list(np.arange(10, 13, 1)) + [fov_to_str(i) for i in np.arange(13, 16, 1)] + [fov_to_str(i) for i in np.arange(17, 26, 1)] + [fov_to_str(i) for i in np.arange(27, 36, 1)] # PC3
fov_lst = [fov_to_str(i) for i in np.arange(1, 2, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 8, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(9, 19, 1)]  # PC9

n_nuclear_convex_dilation = 0
local_size = 150

data = pd.DataFrame(columns=['nuclear', 'FOV',
                             'centroid_nuclear', 'area_nuclear',
                             'radial_curve_nuclear', 'radial_curve_RNAFISH', 'radial_curve_normalized'])

for fov in range(len(fov_lst)):
    print(fov)

    file_name = '%s_%s_Lng_SVCC_Processed001_RAW' % (sample, fov_lst[fov])
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_RNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov_lst[fov]), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))
    else:
        img_nuclear_seg = img_seg

    # measure
    # get local images
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):

        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_hoechst.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_RNAFISH = img_RNAFISH.copy()
        local_RNAFISH = local_RNAFISH[position[0]:position[1], position[2]:position[3]]

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        area_nuclear = local_nuclear_props[0].area

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_RNAFISH = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_RNAFISH, 0.05, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_nuclear, 0.05, 1)

        radial_distribution_relative_r_normalized = \
            list(np.array(radial_distribution_relative_r_RNAFISH) / np.array(
                radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_normalized = dat.nan_replace(radial_distribution_relative_r_normalized)

        data.loc[len(data.index)] = [i, fov,
                                     original_centroid_nuclear, area_nuclear,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_RNAFISH,
                                     radial_distribution_relative_r_normalized]

data.to_csv('%s%s_radial_n%s.txt' % (output_dir, sample, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")