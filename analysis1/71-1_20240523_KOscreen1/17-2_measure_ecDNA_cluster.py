import skimage.io as skio
import napari
import numpy as np
import pandas as pd
import shared.image as ima
import nd2
import shared.dataframe as dat
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B9'
samples = ['B9']
total_fovs = [16]

dshape_factor = 0.145
pixel_size = 300 / 2720  # uM
local_size = 200

data = pd.DataFrame(columns=['sample', 'fov', 'label_nuclear', 'area_nuclear',
                             'n_ecDNA', 'area_ind_ecDNA', 'total_area_ecDNA',
                             'area_ratio_ind_ecDNA', 'total_area_ratio_ecDNA',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half', 'per_AUC'])

pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

pd_thresh = pd.read_csv('%s/%s/17_%s_DNAFISH_thresh.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)
    pd_thresh_s = pd_thresh[pd_thresh['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov + 1, total_fov))
        thresh = pd_thresh_s['DNAFISH_thresh'][fov]
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_hoechst_seg_z = img_hoechst[seg_z]
        img_DNAFISH_seg_z = img_DNAFISH[seg_z]

        img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, s, fov),
                                plugin="tifffile")
        img_seg_ecDNA = skio.imread("%s/%s/17_seg_DNAFISH_tif/%s_%s_seg_DNAFISH.tif" % (output_dir, sample, s, fov),
                                    plugin="tifffile")

        nuclear_props = regionprops(label(img_seg_z), img_seg_z)

        for i in range(len(nuclear_props)):
            print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov, i + 1, len(nuclear_props)))
            original_centroid_nuclear = nuclear_props[i].centroid
            label_nuclear = nuclear_props[i].intensity_mean
            position = ima.img_local_position(img_seg_z, original_centroid_nuclear, local_size)
            local_nuclear_seg = ima.img_local_seg(img_seg_z, position, nuclear_props[i].intensity_mean)
            local_nuclear = img_hoechst_seg_z.copy()[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH = img_DNAFISH_seg_z.copy()[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH_seg = img_seg_ecDNA.copy()[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH_seg[local_nuclear_seg == 0] = 0

            """viewer = napari.Viewer()
            viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 12000])
            viewer.add_image(local_nuclear_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
            napari.run()"""

            # basic measurements
            local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)
            local_DNAFISH_props = regionprops(label(local_nuclear_seg), local_DNAFISH)
            ecDNA_props = regionprops(label(local_DNAFISH_seg), local_DNAFISH)

            area_nuclear = local_nuclear_props[0].area
            # perimeter_nuclear = local_nuclear_props[0].perimeter
            # mean_int_nuclear = local_nuclear_props[0].intensity_mean
            # mean_int_DNAFISH = local_DNAFISH_props[0].intensity_mean
            # circ_nuclear = (4 * math.pi * area_nuclear) / (perimeter_nuclear ** 2)

            if thresh < 1000:
                data.loc[len(data.index)] = [s, fov, label_nuclear, area_nuclear, -1, [-1], -1, [-1], -1, [-1], -1,
                                             [-1], -1, -1]
            else:

                # ecDNA measurements
                n_ecDNA = len(ecDNA_props)
                # centroid_ind_ecDNA = [ecDNA_props[i].centroid for i in range(len(ecDNA_props))]
                area_ind_ecDNA = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
                area_ratio_ind_ecDNA = list(np.array(area_ind_ecDNA) / area_nuclear)
                total_area_ecDNA = sum(area_ind_ecDNA)
                total_area_ratio_ecDNA = sum(area_ratio_ind_ecDNA)

                # mean_int_ind_ecDNA = [ecDNA_props[i].intensity_mean for i in range(len(ecDNA_props))]
                # total_int_ind_ecDNA = list(np.array(area_ind_ecDNA) * np.array(mean_int_ind_ecDNA))
                # total_int_ecDNA = sum(total_int_ind_ecDNA)
                # mean_int_ecDNA = total_int_ecDNA / total_area_ecDNA if total_area_ecDNA != 0 else 0

                # percetage and cumulative curves
                percentage_area_ind_ecDNA = list(np.array(sorted(area_ind_ecDNA, reverse=True)) / total_area_ecDNA)

                cum_percentage_area_ind_ecDNA = dat.list_sum(percentage_area_ind_ecDNA)
                per_AUC = np.sum(cum_percentage_area_ind_ecDNA[:10]) / 10
                cum_percentage_area_n_half = dat.find_pos(0.5, percentage_area_ind_ecDNA)
                cum_area_ind_ecDNA = dat.list_sum(area_ind_ecDNA)
                cum_area_n_half = dat.find_pos(cum_area_ind_ecDNA[-1] / 2, cum_area_ind_ecDNA)

                data.loc[len(data.index)] = [s, fov, label_nuclear, area_nuclear,
                                             n_ecDNA, area_ind_ecDNA, total_area_ecDNA,
                                             area_ratio_ind_ecDNA, total_area_ratio_ecDNA,
                                             cum_percentage_area_ind_ecDNA, cum_percentage_area_n_half,
                                             cum_area_ind_ecDNA, cum_area_n_half, per_AUC]

data['cum_area_ind_ecDNA_filled'] = dat.list_fill_with_last_num(data['cum_area_ind_ecDNA'].tolist())
data.to_csv('%s/%s/17_%s_cluster.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")
