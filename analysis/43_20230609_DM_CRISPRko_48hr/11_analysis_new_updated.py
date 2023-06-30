import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import os
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
sample = 'C2'
total_fov = 16
dshape_factor = 0.0765
n_nuclear_convex_dilation = 4
local_size = 200
rmax = 100
start_fov = 0


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


df_align = pd.read_csv('%s/alignment/%s/%s_alignment_autocheck.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (data_dir, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'FOV', 'label',
                             'n_nuclear_convex_dilation',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear',
                             'mean_int_red', 'mean_int_green',
                             'mean_int_DNAFISH', 'n_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH'])

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s0%s' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_s%s' % (row, sample, fov)
    img_nuclear_bgc = skio.imread("%s/FISH/%s/%s_ch01.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_DNAFISH_bgc = skio.imread("%s/FISH/%s/%s_ch00.tif" % (data_dir, sample, filename), plugin="tifffile")
    img_seg = skio.imread("%s/seg/%s/seg_tif/%s_%s_seg.tif" % (data_dir, sample, sample, fov),
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

        area_nuclear = local_nuclear_props[0].area
        perimeter_nuclear = local_nuclear_props[0].perimeter
        mean_int_nuclear = local_nuclear_props[0].intensity_mean
        mean_int_red = local_red_props[0].intensity_mean
        mean_int_green = local_green_props[0].intensity_mean
        mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
        circ_nuclear = (4 * math.pi * area_nuclear) / (perimeter_nuclear ** 2)

        total_area_ecDNA = []
        total_area_ratio_ecDNA = []
        n_ecDNA = []

        # ecDNA measurements
        for j in range(65):
            thresh = j*1000
            local_DNAFISH_seg = np.zeros_like(local_DNAFISH)
            local_DNAFISH_seg[local_DNAFISH > thresh] = 1
            local_DNAFISH_seg[local_DNAFISH <= thresh+1000] = 0
            # local_DNAFISH_seg = obj.remove_small(local_DNAFISH_seg, 20)
            local_DNAFISH_seg[local_nuclear_seg == 0] = 0
            ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)
            """n_ecDNA_temp = len(ecDNA_props)
            centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
            area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
            area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
            total_area_ecDNA_temp = sum(area_ind_ecDNA)
            total_area_ratio_ecDNA_temp = sum(area_ratio_ind_ecDNA)"""
            total_area_ecDNA_temp = local_DNAFISH_seg.sum()/local_nuclear_seg.sum()
            # total_area_ratio_ecDNA_temp = total_area_ecDNA_temp/area_nuclear

            n_ecDNA.append(0)
            total_area_ecDNA.append(total_area_ecDNA_temp)
            total_area_ratio_ecDNA.append(0 )

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
                                     total_area_ecDNA, total_area_ratio_ecDNA,
                                     radial_distribution_relative_r_nuclear, radial_distribution_relative_r_DNAFISH]

if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
data.to_csv('%s%s/%s_n%s_updated1.txt' % (output_dir, sample, sample, n_nuclear_convex_dilation), index=False, sep='\t')

print("DONE!")