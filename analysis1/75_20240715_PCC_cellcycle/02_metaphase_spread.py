from skimage.filters import threshold_local, sobel
from skimage.morphology import binary_dilation, binary_erosion, dilation
from skimage.measure import label, regionprops
import shared.objects as obj
from skimage.segmentation import watershed
import shared.segmentation as seg
import skimage.io as skio
from skimage.morphology import extrema
import shared.image as ima
from scipy import ndimage
import napari
import shared.display as dis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import nd2

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240715_analysis_PCC_cellcycle/"
data_dir = '%sdata/' % master_folder
save_path = '%sprocessed/' % master_folder

sample = '03_S'
start_fov = 1
total_fov = 120
mode = 'DM_only'  # accept 'call_HSR' or 'DM_only'

###############################################
# ####  PLEASE DO NOT CHANGE BELOW HERE  #### #
###############################################

# PARAMETERS
end_fov = total_fov + 1 - start_fov
if not os.path.exists(save_path):
    os.makedirs(save_path)


# IMAGING ANALYSIS
export_file = '%s%s.txt' % (save_path, sample)
if os.path.exists(export_file):
    data = pd.read_csv(export_file, na_values=['.'], sep='\t')
else:
    data = pd.DataFrame(columns=['nuclear', 'FOV', 'DM_n', 'DM_ind_mean_int', 'DM_ind_area', 'DM_ind_total_int',
                                 'DM_total_int', 'DM_copy', 'HSR_n', 'HSR_ind_mean_int', 'HSR_ind_area',
                                 'HSR_ind_total_int', 'HSR_total_int', 'DM_percentage', 'HSR_percentage'])
# load images
nuclear = 0
for f in range(end_fov):
    fov = f + start_fov
    print(fov)
    if fov < 10:
        img_stack = nd2.imread('%s/%s/%s_00%s.nd2' % (data_dir, sample, sample, fov))
    else:
        img_stack = nd2.imread('%s/%s/%s_0%s.nd2' % (data_dir, sample, sample, fov))
    img_hoechst = img_stack[:, 0, :, :].max(axis=0)
    img_FISH = img_stack[:, 1, :, :].max(axis=0)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_FISH, blending='additive', colormap='green', contrast_limits=[0, img_FISH.max()*0.6])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']

    # analysis on each separate cell
    for i in range(len(shapes.data)):
        nuclear = nuclear + 1
        poly_data = shapes.data[i]

        # ***** generate local images *****

        # reshape sub_masks
        top, left = np.floor(np.min(poly_data, axis=0))
        bottom, right = np.ceil(np.max(poly_data, axis=0))
        top, bottom = np.clip((top, bottom), 0, img_FISH.shape[0] - 1).astype(int)
        left, right = np.clip((left, right), 0, img_FISH.shape[1] - 1).astype(int)
        output_shape = (bottom - top + 1, right - left + 1)
        # generate sub_masks and sub_channels
        sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
        cell_FISH_mask = img_FISH[top:bottom + 1, left:right + 1] * sub_masks
        cell_hoechst_mask = img_hoechst[top:bottom + 1, left:right + 1] * sub_masks
        cell_FISH = img_FISH[top:bottom + 1, left:right + 1]
        cell_hoechst = img_hoechst[top:bottom + 1, left:right + 1]
        # background measurement
        bg_FISH = seg.get_bg_int([cell_FISH])[0]
        bg_hoechst = seg.get_bg_int([cell_hoechst])[0]
        cell_hoechst_bg_corrected = cell_hoechst.astype(float) - np.ones_like(cell_hoechst) * bg_hoechst
        cell_hoechst_bg_corrected[cell_hoechst_bg_corrected < 0] = 0

        # ***** chromosome segmentation *****

        # generate primary chromosome segmentation
        local_hoechst = threshold_local(cell_hoechst, 71)  # original 31 for Aarohi
        chromosome_seg = cell_hoechst_mask > local_hoechst
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_large(chromosome_seg, 10000)  # original 5000 for Aarohi
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_small(chromosome_seg, 100)
        chromosome_seg = ndimage.binary_fill_holes(chromosome_seg)

        # *** filter chromosome segmentation ***
        # measure chromosome mean intensity (hoechst channel)
        chromosome_mean_int_in_hoechst_bg_corrected = np.sum(cell_hoechst_bg_corrected * chromosome_seg) / np.sum(
            chromosome_seg)
        # erose and dilate chromosome segmentation to separate chromosome better
        chromosome_seg_erosion = binary_erosion(chromosome_seg)
        chromosome_seg_erosion = binary_erosion(chromosome_seg_erosion)
        chromosome_seg_label = label(chromosome_seg_erosion)
        for k in range(4):
            chromosome_seg_label = dilation(chromosome_seg_label)
        # remove low intensity ones (presumably ecDNA)
        chromosome_seg_filter = np.zeros_like(cell_hoechst_mask)
        chromosome_props_bg_corrected = regionprops(chromosome_seg_label, cell_hoechst_bg_corrected)
        for j in range(len(chromosome_props_bg_corrected)):
            if chromosome_props_bg_corrected[j].intensity_mean > (chromosome_mean_int_in_hoechst_bg_corrected * 0.5):
                chromosome_seg_filter[chromosome_seg_label == chromosome_props_bg_corrected[j].label] = 1

        # ***** FISH signal segmentation *****

        # measure maximum FISH intensity
        FISH_max = cell_FISH_mask.max()
        # identify local maxima in FISH channel
        maxima_value = 400
        if FISH_max < 4000:
            maxima_value = 250
        maxima = extrema.h_maxima(cell_FISH_mask, maxima_value)
        maxima_outside_chromosome = maxima.copy()
        maxima_outside_chromosome[chromosome_seg_filter == 1] = 0
        elevation_map = sobel(cell_FISH_mask)
        # primary screening for large FISH signal
        markers = np.zeros_like(cell_FISH_mask)
        markers[cell_FISH_mask < (FISH_max / 2)] = 1
        markers[maxima == 1] = 2
        FISH_seg_large = watershed(elevation_map, markers)
        FISH_seg_large_label = obj.label_remove_small(label(FISH_seg_large), 150)
        FISH_seg_large = np.zeros_like(cell_FISH_mask)
        FISH_seg_large[FISH_seg_large_label > 1] = 1
        # adjust FISH_max for images containing large FISH signal
        if FISH_seg_large.max() == 1:
            FISH_max = FISH_max/2
        # identify high intensity regions (mostly on the chromosome)
        markers = np.zeros_like(cell_FISH_mask)
        markers[cell_FISH_mask < (FISH_max / 2)] = 1
        markers[maxima == 1] = 2
        FISH_seg_highInt = watershed(elevation_map, markers)
        FISH_seg_highInt_label = obj.label_resort(obj.label_remove_small(label(FISH_seg_highInt), 2))
        # identify low intensity regions (outside chromosome)
        markers = np.zeros_like(cell_FISH_mask)
        markers[cell_FISH_mask < (FISH_max / 5)] = 1  # 2500
        markers[maxima_outside_chromosome == 1] = 2
        FISH_seg_outsideChromosome = watershed(elevation_map, markers)
        FISH_seg_outsideChromosome_label = obj.label_remove_small(label(FISH_seg_outsideChromosome), 2)
        FISH_seg_outsideChromosome_label = obj.label_remove_large(FISH_seg_outsideChromosome_label, 100)
        # generate total FISH segmentation
        FISH_seg_total = np.zeros_like(cell_FISH_mask)
        FISH_seg_total[FISH_seg_highInt_label > 1] = 1
        FISH_seg_total[FISH_seg_outsideChromosome_label > 1] = 1

        if mode == 'call_HSR':
            # generate HSR segmentation
            # filter for:
            # 1) centroid on the chromosome &
            # 2) area > 70 | neighbour dots (within 12) with dilation of 5 can separate chromosome
            FISH_seg_HSR_label = label(FISH_seg_total.copy())
            HSR_props = regionprops(FISH_seg_HSR_label)
            FISH_seg_HSR = np.zeros_like(cell_FISH_mask)
            for j in range(len(HSR_props)):
                if chromosome_seg_filter[int(HSR_props[j].centroid[0])][int(HSR_props[j].centroid[1])] == 1:
                    if HSR_props[j].area > 70:
                        FISH_seg_HSR[FISH_seg_HSR_label == HSR_props[j].label] = 1
                    else:
                        FISH_seg_small_temp = np.zeros_like(cell_FISH_mask)
                        FISH_seg_small_temp[FISH_seg_HSR_label == HSR_props[j].label] = 1
                        centroid = HSR_props[j].centroid
                        for k in range(len(HSR_props)):
                            centroid_distance = ((HSR_props[k].centroid[0] - centroid[0]) ** 2
                                                 + (HSR_props[k].centroid[1] - centroid[1]) ** 2) ** 0.5
                            if centroid_distance < 12:
                                FISH_seg_small_temp[FISH_seg_HSR_label == HSR_props[k].label] = 1
                        for k in range(5):
                            FISH_seg_small_temp = binary_dilation(FISH_seg_small_temp)
                        chromosome_seg_wo_temp = chromosome_seg_filter.copy()
                        chromosome_seg_wo_temp[FISH_seg_small_temp == 1] = 0
                        if len(regionprops(label(chromosome_seg_filter))) + 1 == \
                                len(regionprops(label(chromosome_seg_wo_temp))):
                            FISH_seg_HSR[FISH_seg_HSR_label == HSR_props[j].label] = 1

            # generate DM segmentation
            FISH_seg_DM = np.zeros_like(cell_FISH_mask)
            FISH_seg_DM = FISH_seg_total.copy()
            FISH_seg_DM[FISH_seg_HSR == 1] = 0

        elif mode == 'DM_only':
            FISH_seg_DM = FISH_seg_total.copy()
            FISH_seg_HSR = np.zeros_like(FISH_seg_DM)

        # manuel correction of DM and HSR
        viewer = napari.Viewer()
        viewer.add_image(cell_hoechst_mask, blending='additive', colormap='blue',
                         contrast_limits=[0, cell_hoechst_mask.max()])
        viewer.add_image(cell_FISH_mask, blending='additive', colormap='green',
                         contrast_limits=[0, cell_FISH_mask.max() * 0.8])
        viewer.add_image(FISH_seg_DM, blending='additive')
        viewer.add_image(FISH_seg_HSR, blending='additive')
        shapes_DM_to_HSR = viewer.add_shapes(name='change DM to HSR', ndim=2)
        shapes_HSR_to_DM = viewer.add_shapes(name='change HSR to DM', ndim=2)
        shapes_DM_remove = viewer.add_shapes(name='DM to be removed', ndim=2)
        shapes_DM_add = viewer.add_shapes(name='DM to be added', ndim=2)
        shapes_HSR_remove = viewer.add_shapes(name='HSR to be removed', ndim=2)
        shapes_HSR_add = viewer.add_shapes(name='HSR to be added', ndim=2)
        napari.run()

        FISH_seg_DM, FISH_seg_HSR = ima.napari_change_between_masks(shapes_DM_to_HSR.data, FISH_seg_DM,
                                                                    FISH_seg_HSR)
        FISH_seg_HSR, FISH_seg_DM = ima.napari_change_between_masks(shapes_HSR_to_DM.data, FISH_seg_HSR,
                                                                    FISH_seg_DM)
        FISH_seg_DM = ima.napari_add_or_remove(shapes_DM_remove.data, 'remove', FISH_seg_DM)
        FISH_seg_DM = ima.napari_add_or_remove(shapes_DM_add.data, 'add', FISH_seg_DM)
        FISH_seg_HSR = ima.napari_add_or_remove(shapes_HSR_remove.data, 'remove', FISH_seg_HSR)
        FISH_seg_HSR = ima.napari_add_or_remove(shapes_HSR_add.data, 'add', FISH_seg_HSR)

        # measure mean chromosome intensity on FISH channel
        chromosome_woFISH_seg = np.zeros_like(cell_FISH_mask)
        chromosome_woFISH_seg = chromosome_seg_filter.copy()
        chromosome_woFISH_seg[FISH_seg_HSR == 1] = 0
        chromosome_bg_in_FISH = np.sum(cell_FISH_mask * chromosome_woFISH_seg) / np.sum(chromosome_woFISH_seg)
        # generate bg corrected FISH image
        cell_FISH_bg_corrected = np.zeros_like(cell_FISH_mask)
        cell_FISH_bg_corrected = cell_FISH_mask.copy().astype(float) \
                                 - (chromosome_seg_filter * (chromosome_bg_in_FISH - bg_FISH)) \
                                 - (np.ones_like(cell_FISH_mask) * bg_FISH)
        cell_FISH_bg_corrected[cell_FISH_bg_corrected < 0] = 0
        cell_FISH_bg_corrected = cell_FISH_bg_corrected.astype(int)

        # measure HSR and DM properties
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
        DM_percentage = DM_total_int / (DM_total_int + HSR_total_int)
        HSR_percentage = HSR_total_int / (DM_total_int + HSR_total_int)

        data.loc[len(data.index)] = [nuclear, fov, len(DM_props), DM_ind_mean_int, DM_ind_area, DM_ind_total_int,
                                     DM_total_int, DM_copy, len(HSR_props), HSR_ind_mean_int, HSR_ind_area,
                                     HSR_ind_total_int, HSR_total_int, DM_percentage, HSR_percentage]
        data.to_csv('%s%s.txt' % (save_path, sample), index=False, sep='\t')

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_hoechst_mask, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst_mask.max()])
        viewer1.add_image(cell_FISH_mask, blending='additive', colormap='green', contrast_limits=[0, cell_FISH_mask.max() * 0.8])
        plt.imsave('%s%s_%s_%s.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.add_image(FISH_seg_HSR, blending='additive')
        plt.imsave('%s%s_%s_%s_HSR.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer2 = napari.Viewer()
        viewer2.add_image(cell_hoechst_mask, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst_mask.max()])
        viewer2.add_image(cell_FISH_mask, blending='additive', colormap='green', contrast_limits=[0, cell_FISH_mask.max() * 0.8])
        viewer2.add_image(FISH_seg_DM, blending='additive')
        plt.imsave('%s%s_%s_%s_DM.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer2))
        viewer2.close()

        viewer3 = napari.Viewer()
        viewer3.add_image(cell_hoechst_mask, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst_mask.max()])
        viewer3.add_image(cell_FISH_mask, blending='additive', colormap='green', contrast_limits=[0, cell_FISH_mask.max() * 0.8])
        viewer3.add_image(chromosome_seg_filter, blending='additive', colormap='yellow')
        plt.imsave('%s%s_%s_%s_chromosome.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer3))
        viewer3.close()

print("DONE!")







