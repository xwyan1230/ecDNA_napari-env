from skimage.filters import threshold_local, sobel
from skimage.morphology import binary_dilation, binary_erosion, dilation
from skimage.measure import label, regionprops
import shared.objects as obj
from skimage.segmentation import watershed
import skimage.io as skio
from skimage.morphology import extrema
import shared.image as ima
import napari
import shared.display as dis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220809_Aarohi_chromosomeRecombination_metaphase" \
                "/08052022_ChrRecomb_KCl_timept/"
sample = 'XY'
save_path = "%sv3_img/%s/" % (master_folder, sample)
if not os.path.exists(save_path):
    os.makedirs(save_path)
start_fov = 6
total_fov = 70 - start_fov

# IMAGING ANALYSIS
export_file = '%s%s.txt' % (save_path, sample)
if os.path.exists(export_file):
    data = pd.read_csv(export_file, na_values=['.'], sep='\t')
else:
    data = pd.DataFrame(columns=['nuclear', 'FOV', 'DM_n', 'DM_ind_mean_int', 'DM_ind_area', 'DM_ind_total_int',
                                 'DM_total_int', 'DM_copy', 'HSR_n', 'HSR_ind_mean_int', 'HSR_ind_area',
                                 'HSR_ind_total_int', 'HSR_total_int', 'HSR_copy', 'DM_percentage', 'HSR_percentage'])
# load images
nuclear = 0
for f in range(total_fov):
    fov = f + start_fov
    print(fov)
    img_hoechst = skio.imread("%s%s/ColoDM_%s_12m_%s_RAW_ch00.tif" % (master_folder, sample, sample, fov), plugin="tifffile")
    img_FISH = skio.imread("%s%s/ColoDM_%s_12m_%s_RAW_ch01.tif" % (master_folder, sample, sample, fov), plugin="tifffile")

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

        # reshape sub_masks
        top, left = np.floor(np.min(poly_data, axis=0))
        bottom, right = np.ceil(np.max(poly_data, axis=0))
        top, bottom = np.clip((top, bottom), 0, img_FISH.shape[0] - 1).astype(int)
        left, right = np.clip((left, right), 0, img_FISH.shape[1] - 1).astype(int)
        output_shape = (bottom - top + 1, right - left + 1)

        # generate sub_masks and sub_channels
        sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
        cell_FISH = img_FISH[top:bottom+1, left:right+1] * sub_masks
        cell_hoechst = img_hoechst[top:bottom+1, left:right+1] * sub_masks

        # generate primary chromosome segmentation
        local_hoechst = threshold_local(cell_hoechst, 31)
        chromosome_seg = cell_hoechst > local_hoechst
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_large(chromosome_seg, 5000)
        chromosome_seg = binary_erosion(chromosome_seg)
        chromosome_seg = obj.remove_small(chromosome_seg, 100)

        # measure chromosome mean intensity
        chromosome_seg_mean_int = np.sum(cell_hoechst * chromosome_seg) / np.sum(chromosome_seg)

        # measure image background intensity of hoechst channel
        img_bg_hoechst = sub_masks.copy()
        img_bg_hoechst[chromosome_seg == 1] = 0
        img_bg_hoechst_int = np.sum(cell_hoechst * img_bg_hoechst) / np.sum(img_bg_hoechst)

        # filter chromosome segmentation
        # erose and dilate chromosome segmentation to separate chromosome better
        chromosome_seg_erosion = binary_erosion(chromosome_seg)
        chromosome_seg_erosion = binary_erosion(chromosome_seg_erosion)
        chromosome_seg_label = label(chromosome_seg_erosion)
        for k in range(4):
            chromosome_seg_label = dilation(chromosome_seg_label)
        # remove low intensity ones (presumably ecDNA)
        chromosome_seg_filter = np.zeros_like(cell_hoechst)
        chromosome_props = regionprops(chromosome_seg_label, cell_hoechst)
        for j in range(len(chromosome_props)):
            if ((chromosome_props[j].intensity_mean - img_bg_hoechst_int) >
                    (chromosome_seg_mean_int - img_bg_hoechst_int) * 0.5):
                chromosome_seg_filter[chromosome_seg_label == chromosome_props[j].label] = 1

        # measure maximum FISH intensity
        FISH_max = cell_FISH.max()
        print(FISH_max)

        # identify local maxima in FISH channel
        maxima_value = 450
        if FISH_max < 4000:
            maxima_value = 250
        maxima = extrema.h_maxima(cell_FISH, maxima_value)
        maxima_outside_chromosome = maxima.copy()
        maxima_outside_chromosome[chromosome_seg_filter == 1] = 0

        # generate elevation map for FISH channel
        elevation_map = sobel(cell_FISH)

        # primary screening for large FISH signal
        markers = np.zeros_like(cell_FISH)
        markers[cell_FISH < (FISH_max / 2)] = 1
        markers[maxima == 1] = 2
        FISH_seg_large = watershed(elevation_map, markers)
        FISH_seg_large_label = obj.label_remove_small(label(FISH_seg_large), 150)
        FISH_seg_large = np.zeros_like(cell_FISH)
        FISH_seg_large[FISH_seg_large_label > 1] = 1

        # adjust FISH_max for images containing large FISH signal
        if FISH_seg_large.max() == 1:
            FISH_max = FISH_max/2
            print(FISH_max)

        # identify high intensity regions (mostly on the chromosome)
        markers = np.zeros_like(cell_FISH)
        markers[cell_FISH < (FISH_max / 2)] = 1
        markers[maxima == 1] = 2
        FISH_seg_highInt = watershed(elevation_map, markers)
        FISH_seg_highInt_label = obj.label_resort(obj.label_remove_small(label(FISH_seg_highInt), 2))

        # identify low intensity regions (outside chromosome)
        markers = np.zeros_like(cell_FISH)
        markers[cell_FISH < (FISH_max / 5)] = 1  # 2500
        markers[maxima_outside_chromosome == 1] = 2
        FISH_seg_outsideChromosome = watershed(elevation_map, markers)
        FISH_seg_outsideChromosome_label = obj.label_remove_small(label(FISH_seg_outsideChromosome), 2)
        FISH_seg_outsideChromosome_label = obj.label_remove_large(FISH_seg_outsideChromosome_label, 100)

        # generate total FISH segmentation
        FISH_seg_total = np.zeros_like(cell_FISH)
        FISH_seg_total[FISH_seg_highInt_label > 1] = 1
        FISH_seg_total[FISH_seg_outsideChromosome_label > 1] = 1

        # generate HSR segmentation
        # filter for:
        # 1) centroid on the chromosome &
        # 2) area > 70 | neighbour dots (within 12) with dilation of 5 can separate chromosome
        FISH_seg_HSR_label = label(FISH_seg_total.copy())
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
        FISH_seg_DM = np.zeros_like(cell_FISH)
        FISH_seg_DM = FISH_seg_total.copy()
        FISH_seg_DM[FISH_seg_HSR == 1] = 0

        # manuel correction of DM and HSR
        viewer = napari.Viewer()
        viewer.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst.max()])
        viewer.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, cell_FISH.max() * 0.8])
        viewer.add_image(FISH_seg_DM, blending='additive')
        viewer.add_image(FISH_seg_HSR, blending='additive')
        shapes_DM_to_HSR = viewer.add_shapes(name='change DM to HSR', ndim=2)
        shapes_HSR_to_DM = viewer.add_shapes(name='change HSR to DM', ndim=2)
        shapes_DM_remove = viewer.add_shapes(name='DM to be removed', ndim=2)
        shapes_DM_add = viewer.add_shapes(name='DM to be added', ndim=2)
        shapes_HSR_remove = viewer.add_shapes(name='HSR to be removed', ndim=2)
        shapes_HSR_add = viewer.add_shapes(name='HSR to be added', ndim=2)
        napari.run()

        FISH_seg_DM, FISH_seg_HSR = ima.napari_change_between_masks(shapes_DM_to_HSR.data, FISH_seg_DM, FISH_seg_HSR)
        FISH_seg_HSR, FISH_seg_DM = ima.napari_change_between_masks(shapes_HSR_to_DM.data, FISH_seg_HSR, FISH_seg_DM)
        FISH_seg_DM = ima.napari_add_or_remove(shapes_DM_remove.data, 'remove', FISH_seg_DM)
        FISH_seg_DM = ima.napari_add_or_remove(shapes_DM_add.data, 'add', FISH_seg_DM)
        FISH_seg_HSR = ima.napari_add_or_remove(shapes_HSR_remove.data, 'remove', FISH_seg_HSR)
        FISH_seg_HSR = ima.napari_add_or_remove(shapes_HSR_add.data, 'add', FISH_seg_HSR)

        # measure image background intensity of FISH channel
        img_bg_FISH = np.zeros_like(cell_FISH)
        img_bg_FISH = sub_masks.copy()
        img_bg_FISH[chromosome_seg_filter == 1] = 0
        img_bg_FISH[FISH_seg_total == 1] = 0
        img_bg_int = np.sum(cell_FISH * img_bg_FISH) / np.sum(img_bg_FISH)

        # measure mean chromosome intensity on FISH channel
        chromosome_woFISH_seg = np.zeros_like(cell_FISH)
        chromosome_woFISH_seg = chromosome_seg_filter.copy()
        chromosome_woFISH_seg[FISH_seg_HSR == 1] = 0
        chromosome_bg_in_FISH = np.sum(cell_FISH * chromosome_woFISH_seg)/np.sum(chromosome_woFISH_seg)

        # generate bg corrected FISH image
        cell_FISH_bg_corrected = np.zeros_like(cell_FISH)
        cell_FISH_bg_corrected = cell_FISH.copy().astype(float) \
                                 - (chromosome_seg_filter * (chromosome_bg_in_FISH - img_bg_int)) \
                                 - (np.ones_like(cell_FISH) * img_bg_int)
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
        if DM_total_int != 0:
            HSR_copy = int(HSR_total_int * DM_copy / DM_total_int)
        else:
            HSR_copy = -1
        DM_percentage = DM_total_int / (DM_total_int + HSR_total_int)
        HSR_percentage = HSR_total_int / (DM_total_int + HSR_total_int)

        data.loc[len(data.index)] = [nuclear, fov, len(DM_props), DM_ind_mean_int, DM_ind_area, DM_ind_total_int,
                                     DM_total_int, DM_copy, len(HSR_props), HSR_ind_mean_int, HSR_ind_area,
                                     HSR_ind_total_int, HSR_total_int, HSR_copy, DM_percentage, HSR_percentage]
        data.to_csv('%s%s.txt' % (save_path, sample), index=False, sep='\t')

        viewer1 = napari.Viewer()
        viewer1.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst.max()])
        viewer1.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, cell_FISH.max() * 0.8])
        plt.imsave('%s%s_%s_%s.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.add_image(FISH_seg_HSR, blending='additive')
        plt.imsave('%s%s_%s_%s_HSR.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer1))
        viewer1.close()

        viewer2 = napari.Viewer()
        viewer2.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst.max()])
        viewer2.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, cell_FISH.max() * 0.8])
        viewer2.add_image(FISH_seg_DM, blending='additive')
        plt.imsave('%s%s_%s_%s_DM.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer2))
        viewer2.close()

        viewer3 = napari.Viewer()
        viewer3.add_image(cell_hoechst, blending='additive', colormap='blue', contrast_limits=[0, cell_hoechst.max()])
        viewer3.add_image(cell_FISH, blending='additive', colormap='green', contrast_limits=[0, cell_FISH.max() * 0.8])
        viewer3.add_image(chromosome_seg_filter, blending='additive', colormap='yellow')
        plt.imsave('%s%s_%s_%s_chromosome.tiff' % (save_path, sample, fov, nuclear), dis.blending(viewer3))
        viewer3.close()

print("DONE!")







