from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis, dilation
from skimage.measure import label, regionprops_table, regionprops
import shared.objects as obj
from skimage.segmentation import watershed
import tifffile as tif
import matplotlib.pyplot as plt
import skimage.io as skio
from skimage.morphology import extrema, remove_small_objects
import napari
import numpy as np
import math
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
sample = 'ATCC'
save_folder = master_folder
total_fov = 76
start_fov = 1

data = pd.DataFrame(columns=['nuclear', 'FOV', 'DM_n', 'DM_ind_mean_int', 'DM_ind_area', 'DM_ind_total_int',
                             'DM_total_int', 'DM_copy', 'HSR_n', 'HSR_ind_mean_int', 'HSR_ind_area',
                             'HSR_ind_total_int', 'HSR_total_int', 'HSR_copy', 'DM_percentage', 'HSR_percentage'])

# IMAGING ANALYSIS
# load images
nuclear = 0
for f in range(total_fov):
    fov = f + start_fov
    img_hoechst = skio.imread("%s%s/ColoDM_%s_12m_%s_RAW_ch00.tif" % (master_folder, sample, sample, fov), plugin="tifffile")
    img_FISH = skio.imread("%s%s/ColoDM_%s_12m_%s_RAW_ch01.tif" % (master_folder, sample, sample, fov), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 20000])
    viewer.add_image(img_FISH, blending='additive', colormap='green', contrast_limits=[0, 10000])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']

    for i in range(len(shapes.data)):
        nuclear = nuclear + 1
        poly_data = shapes.data[i]

        top, left = np.floor(np.min(poly_data, axis=0))
        bottom, right = np.ceil(np.max(poly_data, axis=0))
        top, bottom = np.clip((top, bottom), 0, img_FISH.shape[0] - 1).astype(int)
        left, right = np.clip((left, right), 0, img_FISH.shape[1] - 1).astype(int)
        output_shape = (bottom - top + 1, right - left + 1)

        sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
        cell_FISH = img_FISH[top:bottom+1, left:right+1] * sub_masks
        cell_hoechst = img_hoechst[top:bottom+1, left:right+1] * sub_masks

        local = np.zeros_like(cell_hoechst)
        local = threshold_local(cell_hoechst, 31)  # 21: specific for nucleoli
        chromosome_seg = np.zeros_like(cell_hoechst)
        chromosome_seg = cell_hoechst > local
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_large(chromosome_seg, 5000)
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_small(chromosome_seg, 100)
        chromosome_seg_mean_int = np.sum(cell_hoechst * chromosome_seg) / np.sum(chromosome_seg)

        chromosome_seg_erosion = binary_erosion(chromosome_seg)
        chromosome_seg_erosion = binary_erosion(chromosome_seg_erosion)
        chromosome_seg_label = label(chromosome_seg_erosion)
        for k in range(4):
            chromosome_seg_label = dilation(chromosome_seg_label)

        chromosome_seg_filter = np.zeros_like(cell_hoechst)
        img_bg_hoechst = np.zeros_like(cell_hoechst)
        img_bg_hoechst = sub_masks.copy()
        img_bg_hoechst[chromosome_seg == 1] = 0
        img_bg_hoechst_int = np.sum(cell_hoechst * img_bg_hoechst) / np.sum(img_bg_hoechst)

        chromosome_props = regionprops(chromosome_seg_label, cell_hoechst)

        for j in range(len(chromosome_props)):
            circ = 4 * math.pi * chromosome_props[j].area / (chromosome_props[j].perimeter) ** 2
            if ((chromosome_props[j].intensity_mean - img_bg_hoechst_int) > \
                    (chromosome_seg_mean_int - img_bg_hoechst_int) * 0.5) & (circ < 0.90):
                chromosome_seg_filter[chromosome_seg_label == chromosome_props[j].label] = 1

        FISH_max = cell_FISH.max()
        print(FISH_max)
        print(FISH_max/2)
        print(FISH_max/5)
        maxima = np.zeros_like(cell_FISH)
        maxima = extrema.h_maxima(cell_FISH, 450)  # 700
        maxima_outside_chromosome = np.zeros_like(cell_FISH)
        maxima_outside_chromosome = maxima.copy()
        maxima_outside_chromosome[chromosome_seg_filter == 1] = 0
        elevation_map = np.zeros_like(cell_FISH)
        elevation_map = sobel(cell_FISH)
        markers = np.zeros_like(cell_FISH)
        markers[cell_FISH < (FISH_max/2)] = 1  # 7500
        markers[maxima == 1] = 2
        FISH_seg_highInt = np.zeros_like(cell_FISH)
        FISH_seg_highInt = watershed(elevation_map, markers)
        FISH_seg_highInt_label = np.zeros_like(cell_FISH)
        FISH_seg_highInt_label = obj.label_resort(obj.label_remove_small(label(FISH_seg_highInt), 2))
        FISH_seg_highInt = np.zeros_like(cell_FISH)
        FISH_seg_highInt[FISH_seg_highInt_label > 1] = 1

        markers = np.zeros_like(cell_FISH)
        markers[cell_FISH < (FISH_max/5)] = 1  # 2500
        markers[maxima_outside_chromosome == 1] = 2
        FISH_seg_outsideChromosome = np.zeros_like(cell_FISH)
        FISH_seg_outsideChromosome = watershed(elevation_map, markers)
        FISH_seg_outsideChromosome_label = np.zeros_like(cell_FISH)
        FISH_seg_outsideChromosome_label = obj.label_remove_small(label(FISH_seg_outsideChromosome), 2)
        FISH_seg_outsideChromosome_label = obj.label_remove_large(FISH_seg_outsideChromosome_label, 100)

        FISH_seg_total = np.zeros_like(cell_FISH)
        FISH_seg_total = FISH_seg_highInt.copy()
        FISH_seg_total[FISH_seg_outsideChromosome_label > 1] = 1

        FISH_seg_HSR_label = np.zeros_like(cell_FISH)
        FISH_seg_HSR_label = label(FISH_seg_highInt.copy())
        HSR_props = regionprops(FISH_seg_HSR_label)
        FISH_seg_HSR = np.zeros_like(cell_FISH)
        for j in range(len(HSR_props)):
            if chromosome_seg_filter[int(HSR_props[j].centroid[0])][int(HSR_props[j].centroid[1])] == 1:
                if HSR_props[j].area > 70:
                    FISH_seg_HSR[FISH_seg_HSR_label == HSR_props[j].label] = 1
                else:
                    FISH_seg_small_temp = np.zeros_like(cell_FISH)
                    FISH_seg_small_temp[FISH_seg_HSR_label == HSR_props[j].label] = 1
                    centroid = HSR_props[j].centroid
                    for k in range(len(HSR_props)):
                        centroid_distance = ((HSR_props[k].centroid[0] - centroid[0]) ** 2
                                             + (HSR_props[k].centroid[1] - centroid[1]) ** 2) ** 0.5
                        if centroid_distance < 5:
                            FISH_seg_small_temp[FISH_seg_HSR_label == HSR_props[k].label] = 1
                    for k in range(5):
                        FISH_seg_small_temp = binary_dilation(FISH_seg_small_temp)
                    chromosome_seg_wo_temp = np.zeros_like(cell_FISH)
                    chromosome_seg_wo_temp = chromosome_seg_filter.copy()
                    chromosome_seg_wo_temp[FISH_seg_small_temp == 1] = 0
                    if len(regionprops(label(chromosome_seg_filter))) + 1 == len(regionprops(label(chromosome_seg_wo_temp))):
                        FISH_seg_HSR[FISH_seg_HSR_label == HSR_props[j].label] = 1

        FISH_seg_DM = np.zeros_like(cell_FISH)
        FISH_seg_DM = FISH_seg_total.copy()
        FISH_seg_DM[FISH_seg_HSR == 1] = 0

        img_bg = np.zeros_like(cell_FISH)
        img_bg = sub_masks.copy()
        img_bg[chromosome_seg_filter == 1] = 0
        img_bg[FISH_seg_total == 1] = 0
        img_bg_int = np.sum(cell_FISH * img_bg) / np.sum(img_bg)

        chromosome_woFISH_seg = np.zeros_like(cell_FISH)
        chromosome_woFISH_seg = chromosome_seg_filter.copy()
        chromosome_woFISH_seg[FISH_seg_HSR == 1] = 0
        chromosome_bg_in_FISH = np.sum(cell_FISH * chromosome_woFISH_seg)/np.sum(chromosome_woFISH_seg)

        cell_FISH_bg_corrected = np.zeros_like(cell_FISH)
        cell_FISH_bg_corrected = cell_FISH.copy().astype(float) - (chromosome_seg_filter * (chromosome_bg_in_FISH - img_bg_int)) - \
                                 (np.ones_like(cell_FISH) * img_bg_int)
        cell_FISH_bg_corrected[cell_FISH_bg_corrected < 0] = 0
        cell_FISH_bg_corrected = cell_FISH_bg_corrected.astype(int)

        HSR_props = regionprops(label(FISH_seg_HSR), cell_FISH_bg_corrected)
        HSR_ind_mean_int = [HSR_props[i].intensity_mean for i in range(len(HSR_props))]
        HSR_ind_area = [HSR_props[i].area for i in range(len(HSR_props))]
        HSR_ind_total_int = [HSR_props[i].intensity_mean * HSR_props[i].area for i in range(len(HSR_props))]
        HSR_total_int = sum(HSR_ind_total_int)
        DM_props = regionprops(label(FISH_seg_DM), cell_FISH_bg_corrected)
        DM_ind_mean_int = [DM_props[i].intensity_mean for i in range(len(DM_props))]
        DM_ind_area = [DM_props[i].area for i in range(len(DM_props))]
        DM_ind_total_int = [DM_props[i].intensity_mean * DM_props[i].area for i in range(len(DM_props))]
        DM_total_int = sum(DM_ind_total_int)

        DM_copy = len(DM_props)
        if DM_total_int != 0:
            HSR_copy = int(HSR_total_int * DM_copy / DM_total_int)
        else:
            HSR_copy = -1
        DM_percentage = DM_total_int / (DM_total_int + HSR_total_int)
        HSR_percentage = HSR_total_int / (DM_total_int + HSR_total_int)

        data.loc[len(data.index)] = [nuclear, fov, len(DM_props), DM_ind_mean_int, DM_ind_area, DM_ind_total_int,
                                     DM_total_int, DM_copy, len(HSR_props), HSR_ind_mean_int, HSR_ind_area,
                                     HSR_ind_total_int, HSR_total_int, HSR_copy, DM_percentage, HSR_percentage]

        viewer = napari.Viewer()
        viewer.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 20000])
        viewer.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, 10000])
        viewer.add_image(FISH_seg_HSR)
        viewer.add_image(FISH_seg_DM)
        viewer.add_image(chromosome_seg_filter, blending='additive', colormap='yellow')
        napari.run()


data.to_csv('%s%s.txt' % (save_folder, sample), index=False, sep='\t')

print("DONE!")







