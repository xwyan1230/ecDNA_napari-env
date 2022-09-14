import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis
import shared.objects as obj
import shared.image as ima
from skimage import segmentation
import numpy as np
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import shared.dataframe as dat
import shared.image as ima
import pandas as pd
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220826_neuroblastoma/"
group = 'A'
group_folder = '%s%s/' % (master_folder, group)
save_path = group_folder
sub = 'SZ'
check_lst = ['S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
version = 1

# segmentation
local_factor_nuclear = 151
min_size_nuclear = 3000
max_size_nuclear = 15000
convex_conversion_threshold = 0.85

# other parameters
local_size = 100

data = pd.DataFrame(columns=['sample', 'nuclear',
                             'centroid_nuclear', 'area_nuclear',
                             'n_ecDNA',
                             'total_area_ecDNA', 'total_area_ratio_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA',
                             'max_area_ecDNA', 'max_area_ratio_ecDNA',
                             'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve',
                             'radial_center', 'radial_edge', 'relative_r_area',
                             'angle_curve_nuclear', 'angle_curve_DNAFISH', 'angle_value',
                             'dis_to_hub_area',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'percentage_area_ratio_curve_ecDNA', 'percentage_area_ratio_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half',
                             'cum_area_ratio_ind_ecDNA', 'cum_area_ratio_n_half'])

multi_imgs = [x for x in os.listdir(group_folder)]
if '.DS_Store' in multi_imgs:
    multi_imgs.remove('.DS_Store')

sample_lst = list(set([i.split('.')[0].split('_')[0] for i in multi_imgs]))
sample_lst = [i for i in sample_lst if i[0] in check_lst]

for sample in sample_lst:
    print("Analyzing %s, %s/%s" % (sample, sample_lst.index(sample)+1, len(sample_lst)))

    if os.path.exists("%s%s_b.BMP" % (group_folder, sample)) & os.path.exists("%s%s_r.BMP" % (group_folder, sample)) & os.path.exists("%s%s_g.BMP" % (group_folder, sample)):

        # load images
        img_nuclear = skio.imread("%s%s_b.BMP" % (group_folder, sample))[:, :, 2]
        img_DNAFISH = skio.imread("%s%s_g.BMP" % (group_folder, sample))[:, :, 1]
        img_centromere = skio.imread("%s%s_r.BMP" % (group_folder, sample))[:, :, 0]

        # bg_correction
        bg_val_nuclear = seg.get_bg_int([img_nuclear])[0]
        bg_val_DNAFISH = seg.get_bg_int([img_DNAFISH])[0]
        bg_val_centromere = seg.get_bg_int([img_centromere])[0]
        img_nuclear_bg_corrected = img_nuclear.astype(float) - np.ones_like(img_nuclear) * bg_val_nuclear
        img_nuclear_bg_corrected[img_nuclear_bg_corrected < 0] = 0
        img_DNAFISH_bg_corrected = img_DNAFISH.astype(float) - np.ones_like(img_DNAFISH) * bg_val_DNAFISH
        img_DNAFISH_bg_corrected[img_DNAFISH_bg_corrected < 0] = 0
        img_centromere_bg_corrected = img_centromere.astype(float) - np.ones_like(img_centromere) * bg_val_centromere
        img_centromere_bg_corrected[img_centromere_bg_corrected < 0] = 0

        # nuclear segmentation
        img_nuclear_seg = seg.nuclear_seg(img_nuclear_bg_corrected, local_factor=local_factor_nuclear,
                                          min_size=min_size_nuclear, max_size=max_size_nuclear)
        img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                           threshold=convex_conversion_threshold))
        nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear_bg_corrected)

        img_DNAFISH_seg = np.zeros_like(img_DNAFISH_bg_corrected)
        img_centromere_seg = np.zeros_like(img_centromere_bg_corrected)

        for i in range(len(nuclear_props)):
            original_centroid_nuclear = nuclear_props[i].centroid
            position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
            local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, nuclear_props[i].label)
            local_nuclear = img_nuclear_bg_corrected.copy()
            local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH = img_DNAFISH_bg_corrected.copy()
            local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
            local_centromere = img_centromere_bg_corrected.copy()
            local_centromere = local_centromere[position[0]:position[1], position[2]:position[3]]

            # ecDNA segmentation
            local_DNAFISH_singlet = local_DNAFISH.copy()
            local_DNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
            otsu_threshold_val_local_DNAFISH = threshold_otsu(local_DNAFISH_singlet)

            threshold_min = otsu_threshold_val_local_DNAFISH + (img_DNAFISH_bg_corrected.max() - otsu_threshold_val_local_DNAFISH) / 2
            threshold_min7 = otsu_threshold_val_local_DNAFISH + (img_DNAFISH_bg_corrected.max() - otsu_threshold_val_local_DNAFISH) / 5

            FISH_seg_local = np.zeros_like(local_DNAFISH)

            for k in range(10):
                local = threshold_local(local_DNAFISH, 10*k+7)
                out = local_DNAFISH_singlet > local
                out = binary_erosion(out)
                if k == 0:
                    out = binary_dilation(out)
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH_singlet)
                    for j in range(len(out_props)):
                        temp = np.zeros_like(local_DNAFISH)
                        temp[out_label == out_props[j].label] = 1
                        temp_outer_edge = binary_dilation(temp, disk(4))
                        temp_outer_edge[temp == 1] = 0
                        mean_int_outer_edge = np.sum(local_DNAFISH * temp_outer_edge)/np.sum(temp_outer_edge)
                        if (out_props[j].intensity_mean/mean_int_outer_edge > 1.2) & (out_props[j].area > 5) & \
                                (out_props[j].intensity_mean > threshold_min7):
                            FISH_seg_local[out_label == out_props[j].label] = 1
                else:
                    out_label = label(out)
                    out_props = regionprops(out_label, local_DNAFISH_singlet)
                    for j in range(len(out_props)):
                        if (out_props[j].intensity_mean > threshold_min) & (out_props[j].area > 5):
                            FISH_seg_local[out_label == out_props[j].label] = 1

            bg_val = otsu_threshold_val_local_DNAFISH * 3
            extreme_val = local_DNAFISH_singlet.max() * 2 / otsu_threshold_val_local_DNAFISH
            maxima = extrema.h_maxima(local_DNAFISH, extreme_val)
            elevation_map = sobel(local_DNAFISH)
            markers = np.zeros_like(local_DNAFISH)
            markers[local_DNAFISH_singlet < bg_val] = 1
            markers[maxima == 1] = 2
            seg_wat = segmentation.watershed(elevation_map, markers)
            seg_wat_filter = obj.label_remove_small(label(seg_wat), min_size=6)
            FISH_seg_watershed = np.zeros_like(seg_wat)
            FISH_seg_watershed[seg_wat_filter > 1] = 1

            FISH_seg = FISH_seg_watershed.copy()
            FISH_seg[FISH_seg_local == 1] = 1
            FISH_seg[local_nuclear_seg_convex == 0] = 0

            # centromere segmentation
            local_centromere_singlet = local_centromere.copy()
            local_centromere_singlet[local_nuclear_seg_convex == 0] = 0
            otsu_threshold_val_local_centromere = threshold_otsu(local_centromere_singlet)

            threshold_min = otsu_threshold_val_local_centromere + (
                        img_centromere_bg_corrected.max() - otsu_threshold_val_local_centromere) / 2
            threshold_min7 = otsu_threshold_val_local_centromere + (
                        img_centromere_bg_corrected.max() - otsu_threshold_val_local_centromere) / 5

            centromere_seg_local = np.zeros_like(local_centromere)

            for k in range(10):
                local = threshold_local(local_centromere, 10 * k + 7)
                out = local_centromere_singlet > local
                out = binary_erosion(out)
                if k == 0:
                    out = binary_dilation(out)
                    out_label = label(out)
                    out_props = regionprops(out_label, local_centromere_singlet)
                    for j in range(len(out_props)):
                        temp = np.zeros_like(local_centromere)
                        temp[out_label == out_props[j].label] = 1
                        temp_outer_edge = binary_dilation(temp, disk(4))
                        temp_outer_edge[temp == 1] = 0
                        mean_int_outer_edge = np.sum(local_centromere * temp_outer_edge) / np.sum(temp_outer_edge)
                        if (out_props[j].intensity_mean / mean_int_outer_edge > 1.2) & (out_props[j].area > 5) & \
                                (out_props[j].intensity_mean > threshold_min7):
                            centromere_seg_local[out_label == out_props[j].label] = 1
                else:
                    out_label = label(out)
                    out_props = regionprops(out_label, local_centromere_singlet)
                    for j in range(len(out_props)):
                        if (out_props[j].intensity_mean > threshold_min) & (out_props[j].area > 5):
                            centromere_seg_local[out_label == out_props[j].label] = 1

            bg_val = otsu_threshold_val_local_centromere * 3
            extreme_val = local_centromere_singlet.max() * 2 / otsu_threshold_val_local_centromere
            maxima = extrema.h_maxima(local_centromere, extreme_val)
            elevation_map = sobel(local_centromere)
            markers = np.zeros_like(local_centromere)
            markers[local_centromere_singlet < bg_val] = 1
            markers[maxima == 1] = 2
            seg_wat = segmentation.watershed(elevation_map, markers)
            seg_wat_filter = obj.label_remove_small(label(seg_wat), min_size=6)
            centromere_seg_watershed = np.zeros_like(seg_wat)
            centromere_seg_watershed[seg_wat_filter > 1] = 1

            centromere_seg = centromere_seg_watershed.copy()
            centromere_seg[centromere_seg_local == 1] = 1
            centromere_seg[local_nuclear_seg_convex == 0] = 0

            if (np.sum(FISH_seg)/np.sum(local_nuclear_seg_convex) <= 0.95) & (np.sum(centromere_seg)/np.sum(local_nuclear_seg_convex) <= 0.95):
                img_DNAFISH_seg = ima.image_paste_to(img_DNAFISH_seg, FISH_seg, [int(original_centroid_nuclear[0]-100), int(original_centroid_nuclear[1]-100)])
                img_centromere_seg = ima.image_paste_to(img_centromere_seg, centromere_seg, [int(original_centroid_nuclear[0]-100), int(original_centroid_nuclear[1]-100)])

                # basic measurements
                local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
                ecDNA_props = regionprops(label(FISH_seg))

                area_nuclear = local_nuclear_props[0].area

                # ecDNA measurements
                n_ecDNA = len(ecDNA_props)
                centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
                area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
                area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
                total_area_ecDNA = sum(area_ind_ecDNA)
                total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)
                max_area_ecDNA = max(area_ind_ecDNA + [0])
                max_area_ratio_ecDNA = max(area_ratio_ind_ecDNA + [0])

                # percentage and cumulative curves
                percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True)) / total_area_ecDNA)
                percentage_area_ratio_ind_ecDNA = list(np.array(sorted(area_ratio_ind_ecDNA, reverse=True)) /
                                                       total_area_ratio_ecDNA)

                cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
                cum_percentage_area_n_half = dat.find_pos(0.5, percentage_area_ind_ecDNA)
                cum_percentage_area_ratio_ind_ecDNA = dat.list_sum(percentage_area_ratio_ind_ecDNA)
                cum_percentage_area_ratio_n_half = dat.find_pos(0.5, percentage_area_ratio_ind_ecDNA)

                cum_area_ind_ecDNA = dat.list_sum(area_ind_ecDNA)
                cum_area_n_half = dat.find_pos(cum_area_ind_ecDNA[-1] / 2, cum_area_ind_ecDNA)
                cum_area_ratio_ind_ecDNA = dat.list_sum(area_ratio_ind_ecDNA)
                cum_area_ratio_n_half = dat.find_pos(cum_area_ratio_ind_ecDNA[-1] / 2, cum_area_ratio_ind_ecDNA)

                # distance from the hub
                dis_to_hub_area = 0

                if n_ecDNA == 0:
                    dis_to_hub_area = -1
                elif n_ecDNA > 1:
                    ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'centroid': centroid_ind_ecDNA})

                    ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
                                                                       ignore_index=True)
                    ind_ecDNA_sort_area['dis'] = \
                        [((ind_ecDNA_sort_area['centroid'][i][0] - ind_ecDNA_sort_area['centroid'][0][0]) ** 2 +
                          (ind_ecDNA_sort_area['centroid'][i][1] - ind_ecDNA_sort_area['centroid'][0][1]) ** 2) ** 0.5
                         for i in range(len(ind_ecDNA_sort_area))]

                    for ind in range(n_ecDNA - 1):
                        dis_to_hub_area = dis_to_hub_area + ind_ecDNA_sort_area['area'][ind + 1] / total_area_ecDNA * \
                                          ind_ecDNA_sort_area['dis'][ind + 1]

                # radial distribution
                local_nuclear_centroid = local_nuclear_props[0].centroid
                _, local_edge_distance_map = medial_axis(local_nuclear_seg_convex, return_distance=True)
                local_centroid_distance_map = ima.distance_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
                local_centroid_distance_map[local_nuclear_seg_convex == 0] = 0
                local_edge_distance_map[local_nuclear_seg_convex == 0] = -1
                local_relative_r_map = local_centroid_distance_map / (local_centroid_distance_map + local_edge_distance_map)

                radial_distribution_relative_r_DNAFISH = \
                    ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map, local_DNAFISH, 0.01,
                                                              1)
                radial_distribution_relative_r_nuclear = \
                    ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map, local_nuclear, 0.01,
                                                              1)

                radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
                radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

                radial_curve = list(np.array(radial_distribution_relative_r_DNAFISH_smooth)/np.array(radial_distribution_relative_r_nuclear_smooth))

                radial_subtract = np.array(radial_distribution_relative_r_DNAFISH_smooth) - \
                                  np.array(radial_distribution_relative_r_nuclear_smooth)
                radial_subtract_center = np.mean(radial_subtract[0:40])
                radial_subtract_edge = np.mean(radial_subtract[40:80])

                # relative_r
                relative_r_area = 0

                if n_ecDNA == 0:
                    relative_r_area = -1
                elif n_ecDNA >= 1:
                    ind_ecDNA = pd.DataFrame({'area': area_ind_ecDNA, 'centroid': centroid_ind_ecDNA})

                    ind_ecDNA_sort_area = ind_ecDNA.copy().sort_values(by='area', axis=0, ascending=False, inplace=False,
                                                                       ignore_index=True)
                    ind_ecDNA_sort_area['dis'] = \
                        [local_relative_r_map[int(ind_ecDNA_sort_area['centroid'][i][0])][
                             int(ind_ecDNA_sort_area['centroid'][i][1])]
                         for i in range(len(ind_ecDNA_sort_area))]

                    for ind in range(n_ecDNA):
                        relative_r_area = relative_r_area + ind_ecDNA_sort_area['area'][ind] / total_area_ecDNA * \
                                          ind_ecDNA_sort_area['dis'][ind]

                # angle distribution
                local_angle_map = ima.angle_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
                angle_distribution_DNAFISH = \
                    ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map, local_DNAFISH, 1, 360)
                angle_distribution_nuclear = \
                    ima.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map, local_nuclear, 1, 360)

                angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
                angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
                angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
                    dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)
                angle_value = angle_distribution_DNAFISH_smooth_centered[179]

                data.loc[len(data.index)] = [sample, i, original_centroid_nuclear, area_nuclear, n_ecDNA,
                                             total_area_ecDNA, total_area_ratio_ecDNA, area_ind_ecDNA, area_ratio_ind_ecDNA,
                                             max_area_ecDNA, max_area_ratio_ecDNA,
                                             radial_distribution_relative_r_nuclear_smooth,
                                             radial_distribution_relative_r_DNAFISH_smooth, radial_curve,
                                             radial_subtract_center, radial_subtract_edge, relative_r_area,
                                             angle_distribution_nuclear_smooth_centered,
                                             angle_distribution_DNAFISH_smooth_centered, angle_value,
                                             dis_to_hub_area,
                                             cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                             cum_percentage_area_ratio_ind_ecDNA, cum_percentage_area_ratio_n_half,
                                             cum_area_ind_ecDNA, cum_area_n_half,
                                             cum_area_ratio_ind_ecDNA, cum_area_ratio_n_half]

            """viewer = napari.Viewer()
            viewer.add_image(local_nuclear, blending='additive', colormap='blue')
            viewer.add_image(local_DNAFISH, blending='additive', colormap='green',
                             contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
            viewer.add_image(FISH_seg_local, blending='additive')
            viewer.add_image(FISH_seg, blending='additive')
            napari.run()"""

        # viewer
        """viewer = napari.Viewer()
        viewer.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
        viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
        # viewer.add_image(img_centromere_bg_corrected, blending='additive', colormap='red')
        viewer.add_image(img_nuclear_seg_convex, blending='additive')
        plt.imsave('%s%s_nuclei.tiff' % (save_path, sample), dis.blending(viewer))
        viewer.close()
    
        viewer1 = napari.Viewer()
        viewer1.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
        viewer1.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
        viewer1.add_image(img_DNAFISH_seg, blending='additive')
        plt.imsave('%s%s_DNAFISH.tiff' % (save_path, sample), dis.blending(viewer1))
        viewer1.close()"""

        viewer2 = napari.Viewer()
        viewer2.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
        viewer2.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
        viewer2.add_image(img_centromere_bg_corrected, blending='additive', colormap='red')
        plt.imsave('%s%s_img.tiff' % (save_path, sample), dis.blending(viewer2))
        viewer2.close()

        viewer3 = napari.Viewer()
        viewer3.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue')
        viewer3.add_image(img_DNAFISH_seg, blending='additive', colormap='green')
        viewer3.add_image(img_centromere_seg, blending='additive', colormap='red')
        plt.imsave('%s%s_seg.tiff' % (save_path, sample), dis.blending(viewer3))
        viewer3.close()

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'])
data['cum_area_ratio_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ratio_ind_ecDNA'])

data.to_csv('%s%s_v%s_%s.txt' % (master_folder, group, version, sub), index=False, sep='\t')

print("DONE!")
