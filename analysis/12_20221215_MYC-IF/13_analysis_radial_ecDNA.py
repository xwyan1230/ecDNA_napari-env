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

sample = 'MYC-new'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_MYC_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 1
local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV',
                             'centroid_nuclear', 'area_nuclear', 'mean_int_nuclear', 'MYC_mean',
                             'mean_int_DNAFISH', 'n_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized',
                             'radial_curve_MYC'])

for fov in range(img_MYC_stack.shape[0]):
    print(fov)
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_MYC = img_MYC_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    img_nuclear = np.concatenate([np.zeros(shape=[3, 3144], dtype=np.uint16), img_nuclear], axis=0)[:3144, :3144]
    img_DNAFISH = np.concatenate([np.zeros(shape=[6, 3144], dtype=np.uint16), img_DNAFISH], axis=0)[:3144, :3144]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_seg = np.concatenate([np.zeros(shape=[3, 3144], dtype=np.uint16), img_seg], axis=0)[:3144, :3144]
    img_ecSeg = skio.imread("%s%s/seg_tif/%s_%s_ecseg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_ecSeg = np.concatenate([np.zeros(shape=[6, 3144], dtype=np.uint16), img_ecSeg], axis=0)[:3144, :3144]

    """viewer = napari.Viewer()
    viewer.add_image(img_nuclear_bgc, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bgc.max()])
    viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
    napari.run()"""

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
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_MYC = img_MYC.copy()
        local_MYC = local_MYC[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_nuclear_seg_MYC = ima.img_local_seg(img_seg, position, nuclear_props[i].label)

        """viewer = napari.Viewer()
        viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, local_nuclear.max()])
        viewer.add_image(local_mCherry, blending='additive', colormap='red', contrast_limits=[0, local_mCherry.max()])
        viewer.add_image(local_nuclear_seg_mCherry, blending='additive', colormap='green', contrast_limits=[0, 1])
        viewer.add_image(local_nuclear_seg, blending='additive', contrast_limits=[0,1])
        napari.run()"""

        # basic measurements
        local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
        local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
        local_IF_props = regionprops(label(local_nuclear_seg_MYC), local_MYC)
        ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

        area_nuclear = local_nuclear_props[0].area
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        mean_int_IF = local_IF_props[0].intensity_mean
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean

        # ecDNA measurements
        n_ecDNA = len(ecDNA_props)
        centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
        area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
        area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
        total_area_ecDNA = sum(area_ind_ecDNA)
        total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

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
        radial_distribution_relative_r_MYC = \
            ima.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map,
                                                      local_MYC, 0.05, 1)

        radial_distribution_relative_r_normalized = \
            list(np.array(radial_distribution_relative_r_DNAFISH) / np.array(
                radial_distribution_relative_r_nuclear))
        radial_distribution_relative_r_normalized = dat.nan_replace(radial_distribution_relative_r_normalized)

        data.loc[len(data.index)] = [i, fov,
                                     original_centroid_nuclear, area_nuclear, mean_int_nuclear, mean_int_IF,
                                     mean_int_DNAFISH, n_ecDNA,
                                     total_area_ecDNA, total_area_ratio_ecDNA,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH,
                                     radial_distribution_relative_r_normalized, radial_distribution_relative_r_MYC]

data.to_csv('%s%s_radial_n%s.txt' % (output_dir, sample, n_nuclear_convex_dilation), index=False, sep='\t')
print("DONE!")