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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-old'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_MYC_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 2
local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV', 'centroid_nuclear',
                             'radial_curve_nuclear', 'radial_curve_MYC', 'radial_curve_normalized'])

for fov in range(img_MYC_stack.shape[0]):
    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_MYC = img_MYC_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    # measure
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
        local_mCherry = img_MYC.copy()
        local_mCherry = local_mCherry[position[0]:position[1], position[2]:position[3]]
        local_nuclear_seg_mCherry = ima.img_local_seg(img_seg, position, nuclear_props[i].label)

        local_nuclear_props = regionprops(label(local_nuclear_seg_mCherry), local_nuclear)

        # radial distribution
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg_mCherry, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg_mCherry,
                                                                  local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg_mCherry == 0] = 0
        local_edge_distance_map[local_nuclear_seg_mCherry == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)

        radial_distribution_relative_r_MYC = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg_mCherry, local_relative_r_map,
                                                      local_mCherry,
                                                      0.1, 1)
        radial_distribution_relative_r_nuclear = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg_mCherry, local_relative_r_map,
                                                      local_nuclear,
                                                      0.1, 1)

        radial_distribution_relative_r_normalized = \
            list(np.array(radial_distribution_relative_r_MYC) / np.array(
                radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_normalized = dat.nan_replace(radial_distribution_relative_r_normalized)

        data.loc[len(data.index)] = [i, fov, original_centroid_nuclear,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_MYC,
                                     radial_distribution_relative_r_normalized]

data.to_csv('%s%s_MYC.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")