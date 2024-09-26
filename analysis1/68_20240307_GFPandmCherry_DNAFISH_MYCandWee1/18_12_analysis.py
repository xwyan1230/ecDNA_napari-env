import skimage.io as skio
import pandas as pd
from skimage.measure import label, regionprops
import shared.image as ima
from skimage.morphology import medial_axis
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import imutils
import cv2
import math
import os
import napari

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_fov_pixel = 180
# 766nm for Keyence 10x
# 360nm for 512x512 at 63x
# 58.6nm for 3144x3144 at 63x (0.0765)
# 180nm for 1024x1024 at 63x (0.2349)
# 720nm for 256x256 at 63x (0.9399)

# NO NEED TO CHANGE
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%sfigures/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s' % (sample, treatment, sample)
# filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s' % sample

# DO NOT CHANGE


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


total_fov = 16
start_fov = 0
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
local_size = 100

topleft_final = pd.read_csv('%s/%s/align_topleft_final_refine_%s.txt' % (align_dir, sample, sample), na_values=['.'], sep='\t')
topleft_final['topleft_target'] = [dat.str_to_float(x) for x in topleft_final['topleft_target']]

img_before_GFP = skio.imread("%s/%s/%s_beforeFISH_GFP_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/%s/%s_beforeFISH_mCherry_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'FOV',
                             'centroid_nuclear', 'area_nuclear', 'circ_nuclear', 'mean_int_nuclear',
                             'mean_int_red', 'mean_int_green',
                             'mean_int_DNAFISH', 'n_ecDNA', 'total_area_ecDNA', 'total_area_ratio_ecDNA',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r'])

for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    img_nuclear_ori = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov + 1),
                                  plugin="tifffile")
    img_DNAFISH_ori = skio.imread("%s/DNAFISH/%s/%s_%s_ch00.tif" % (data_dir, sample, filename, fov + 1),
                                  plugin="tifffile")
    img_nuclear = cv2.flip(imutils.rotate(img_nuclear_ori, angle=-90), 0)
    img_DNAFISH = cv2.flip(imutils.rotate(img_DNAFISH_ori, angle=-90), 0)
    img_seg = skio.imread("%s/%s/seg_tif/%s_%s_seg.tif" % (data_dir1, sample, sample, fov + 1), plugin="tifffile")
    img_ecSeg = skio.imread("%s/%s/seg_tif/%s_%s_ecseg.tif" % (data_dir1, sample, sample, fov+1), plugin="tifffile")

    s = img_nuclear.shape[0] * dshape_factor
    topleft = topleft_final['topleft_target'][fov]
    img_before_GFP_cut = img_before_GFP.copy()[int(topleft[1]):int(topleft[1] + s) + 5, int(topleft[0]):int(topleft[0] + s) + 5]
    s1 = img_before_GFP_cut.shape[0]
    img_before_GFP_cut_resize = cv2.resize(img_before_GFP_cut, dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                           interpolation=cv2.INTER_AREA)[:img_nuclear.shape[0], :img_nuclear.shape[1]]
    img_before_mCherry_cut = img_before_mCherry.copy()[int(topleft[1]):int(topleft[1] + s) + 5, int(topleft[0]):int(topleft[0] + s) + 5]
    img_before_mCherry_cut_resize = cv2.resize(img_before_mCherry_cut, dsize=(int(s1 * 1 / dshape_factor), int(s1 * 1 / dshape_factor)),
                                               interpolation=cv2.INTER_AREA)[:img_nuclear.shape[0], :img_nuclear.shape[1]]
    viewer = napari.Viewer()
    viewer.add_image(img_before_mCherry_cut_resize, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_before_GFP_cut_resize, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    plt.imsave('%s/%s/%s_%s.tiff' % (align_dir, sample, sample, fov), dis.blending(viewer))
    # viewer.add_image(img_seg, blending='additive', colormap='green', contrast_limits=[0, 20])
    viewer.close()

    nuclear_props = regionprops(img_seg)

    for i in range(len(nuclear_props)):
        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_red = img_before_mCherry_cut_resize.copy()[position[0]:position[1], position[2]:position[3]]
        local_green = img_before_GFP_cut_resize.copy()[position[0]:position[1], position[2]:position[3]]

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
        # centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
        area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
        area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
        total_area_ecDNA = sum(area_ind_ecDNA)
        total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

        # radial measurements
        local_nuclear_centroid = local_nuclear_props[0].centroid
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map / (
                local_centroid_distance_map + local_edge_distance_map)
        DNAFISH_seg_label = img_to_pixel_int(local_nuclear_seg, local_DNAFISH_seg)
        int_r_to_edge = img_to_pixel_int(local_nuclear_seg, local_edge_distance_map)
        int_relative_r = img_to_pixel_int(local_nuclear_seg, local_relative_r_map)

        data.loc[len(data.index)] = [i, fov+1, original_centroid_nuclear, area_nuclear, circ_nuclear, mean_int_nuclear,
                                     mean_int_red, mean_int_green, mean_int_DNAFISH, n_ecDNA, total_area_ecDNA,
                                     total_area_ratio_ecDNA, DNAFISH_seg_label, int_r_to_edge, int_relative_r]

if not os.path.exists("%s/" % output_dir):
    os.makedirs("%s/" % output_dir)
data.to_csv('%s%s.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")