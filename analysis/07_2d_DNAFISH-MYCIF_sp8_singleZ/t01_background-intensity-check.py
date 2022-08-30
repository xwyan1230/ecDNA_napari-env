from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, binary_dilation, disk, binary_erosion
from skimage.filters import threshold_otsu, try_all_threshold
import pandas as pd
import numpy as np
import shared.image as img
import skimage.io as skio
import shared.image as ima
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.segmentation as seg
import shared.math as mat
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220822_sp8_POLR3Dko_JunDMandHSR_singleZ/"
sample = 'HSR_Jun_POLR3Dko'
save_path = master_folder
start_fov = 0
total_fov = 30
version = 1

local_size = 200
rmax = 100

data = pd.DataFrame(columns=['nuclear', 'FOV',
                             'centroid_nuclear', 'area_nuclear',
                             'bg_int_nuclear', 'mean_int_nuclear', 'total_int_nuclear',
                             'bg_int_IF', 'mean_int_IF', 'total_int_IF',
                             'bg_int_DNAFISH', 'mean_int_DNAFISH', 'total_int_DNAFISH',
                             'n_ecDNA',
                             'mean_int_ecDNA', 'total_int_ecDNA', 'mean_int_ind_ecDNA', 'total_int_ind_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                             'max_area_ecDNA', 'max_area_ratio_ecDNA',
                             'g', 'g_value',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_center', 'radial_edge',
                             'relative_r_area', 'relative_r_int',
                             'angle_curve_nuclear', 'angle_curve_DNAFISH', 'angle_value',
                             'dis_to_hub_area', 'dis_to_hub_int',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'percentage_int_curve_ecDNA', 'percentage_int_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half',
                             'cum_int_ind_ecDNA', 'cum_int_n_half'])

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, fov %s/%s" % (sample, fov+1, total_fov))
    if fov < 10:
        file_prefix = '20220822_sp8_POLR3Dko_JunDMandHSR_singleZ_%s_s0%s' % (sample, fov)
    else:
        file_prefix = '20220822_sp8_POLR3Dko_JunDMandHSR_singleZ_%s_s%s' % (sample, fov)

    # load images
    img_nuclear = skio.imread("%s%s/%s_ch00.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s_ch02.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_IF = skio.imread("%s%s/%s_ch01.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_nuclear_seg = skio.imread("%s%s/%s_seg.tif" % (master_folder, sample, file_prefix), plugin="tifffile")
    img_DNAFISH_seg = skio.imread("%s%s/%s_ecSeg.tif" % (master_folder, sample, file_prefix), plugin="tifffile")

    # background measurement and background correction
    bg_int_nuclear = seg.get_bg_int([img_nuclear])[0]
    bg_int_DNAFISH = seg.get_bg_int([img_DNAFISH])[0]
    bg_int_IF = seg.get_bg_int([img_IF])[0]
    bg_otsu_nuclear = threshold_otsu(img_nuclear)
    bg_otsu_DNAFISH = threshold_otsu(img_DNAFISH)
    bg_otsu_IF = threshold_otsu(img_IF)
    img_nuclear_bgc = img_nuclear > bg_int_nuclear
    img_DNAFISH_bgc = img_DNAFISH > bg_int_DNAFISH
    img_IF_bgc = img_IF > bg_int_IF
    img_nuclear_otsu = img_nuclear > bg_otsu_nuclear
    img_DNAFISH_otsu = img_DNAFISH > bg_otsu_DNAFISH
    img_IF_otsu = img_IF > bg_otsu_IF

    img_nuclear_otsu = binary_dilation(img_nuclear_otsu)
    img_nuclear_otsu = binary_erosion(img_nuclear_otsu, disk(2))
    img_nuclear_otsu = binary_dilation(img_nuclear_otsu)
    img_DNAFISH_otsu = binary_dilation(img_DNAFISH_otsu)
    img_DNAFISH_otsu = binary_erosion(img_DNAFISH_otsu, disk(2))
    img_DNAFISH_otsu = binary_dilation(img_DNAFISH_otsu)
    img_IF_otsu = binary_dilation(img_IF_otsu)
    img_IF_otsu = binary_erosion(img_IF_otsu, disk(2))
    img_IF_otsu = binary_dilation(img_IF_otsu)

    bg = img_nuclear_otsu.copy()
    bg[img_DNAFISH_otsu == 1] = 1
    bg[img_IF_otsu == 1] = 1
    bg = binary_dilation(bg, disk(50))

    background = np.ones_like(bg)
    background[bg == 1] = 0

    bg_nuclear = np.sum(background * img_nuclear)/np.sum(background)
    bg_DNAFISH = np.sum(background * img_DNAFISH)/np.sum(background)
    bg_IF = np.sum(background * img_IF)/np.sum(background)

    print(bg_int_nuclear)
    print(bg_int_DNAFISH)
    print(bg_int_IF)
    print(bg_otsu_nuclear)
    print(bg_otsu_DNAFISH)
    print(bg_otsu_IF)
    print(bg_nuclear)
    print(bg_DNAFISH)
    print(bg_IF)

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear_seg)
    # viewer.add_image(img_nuclear_seg_dilation)
    # viewer.add_image(img_nuclear_bgc)
    # viewer.add_image(img_DNAFISH_bgc)
    # viewer.add_image(img_IF_bgc)
    viewer.add_image(img_nuclear_otsu)
    viewer.add_image(img_DNAFISH_otsu)
    viewer.add_image(img_IF_otsu)
    viewer.add_image(bg)
    napari.run()