import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import tifffile as tif
import shared.image as ima
import shared.dataframe as dat
from skimage.morphology import medial_axis
import nd2
import os
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'B4'
samples = ['B4_1_9pos', 'B4_2_7pos']
total_fovs = [9, 7]

### 14
print("Running 14...")
dshape_factor = 0.145
pixel_size = 300/2720  # uM

img_before_GFP = imutils.rotate(skio.imread("%s/%s/03_%s_GFP_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)
img_before_mCherry = imutils.rotate(skio.imread("%s/%s/03_%s_mCherry_final.tif" % (output_dir, sample, sample), plugin="tifffile"), angle=180)

align = pd.read_csv('%s/%s/08_%s_alignment.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

data = pd.DataFrame()

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    align_s = align[align['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)
        img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")

        dsize = int(img_hoechst_merge.shape[1] * dshape_factor)
        x = int(align_s['topleft_x'][fov])
        y = int(align_s['topleft_y'][fov])
        img_before_GFP_fov = img_before_GFP[y:(y+dsize), x:(x+dsize)]
        img_before_mCherry_fov = img_before_mCherry[y:(y+dsize), x:(x+dsize)]
        img_before_GFP_fov_resize = cv2.resize(img_before_GFP_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)
        img_before_mCherry_fov_resize = cv2.resize(img_before_mCherry_fov, dsize=(img_seg_total.shape[0], img_seg_total.shape[1]), interpolation=cv2.INTER_AREA)

        label_props = regionprops(label(img_seg_total), img_seg_total)
        GFP_props = regionprops(label(img_seg_total), img_before_GFP_fov_resize)
        mCherry_props = regionprops(label(img_seg_total), img_before_mCherry_fov_resize)

        data_temp = pd.DataFrame()
        data_temp['sample'] = [s] * len(label_props)
        data_temp['fov'] = [fov] * len(label_props)
        data_temp['label_mean_int'] = [label_props[i].intensity_mean for i in range(len(label_props))]
        data_temp['GFP'] = [GFP_props[i].intensity_mean for i in range(len(label_props))]
        data_temp['mCherry'] = [mCherry_props[i].intensity_mean for i in range(len(label_props))]
        data = pd.concat([data, data_temp], axis=0)

        if not os.path.exists("%s/%s/14_red_green_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/14_red_green_color/" % (output_dir, sample))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
        viewer.add_image(img_before_GFP_fov_resize, blending='additive', colormap='green', contrast_limits=[0, 40000])
        viewer.add_image(img_before_mCherry_fov_resize, blending='additive', colormap='red', contrast_limits=[0, 40000])
        plt.imsave("%s%s/14_red_green_color/%s_%s_red_green_color.tiff" % (output_dir, sample, s, fov),
                   dis.blending(viewer))
        viewer.close()

data.to_csv('%s/%s/14_%s_red_green.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 16
print("Running 16...")
data = pd.DataFrame()
df_bg = pd.read_csv('%s/bg.txt' % output_dir, na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        img_DNAFISH_sum = img_DNAFISH.astype(int).sum(axis=0)

        img_seg_total = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        label_props = regionprops(label(img_seg_total), img_seg_total)
        DNAFISH_props = regionprops(label(img_seg_total), img_DNAFISH_sum)
        bg_merge = df_bg[df_bg['sample'] == s]['bg_merge'].tolist()[0]

        data_temp = pd.DataFrame()
        data_temp['sample'] = [s] * len(label_props)
        data_temp['fov'] = [fov] * len(label_props)
        data_temp['label_mean_int'] = [label_props[i].intensity_mean for i in range(len(label_props))]
        data_temp['nuclear_area_merge'] = [label_props[i].area for i in range(len(label_props))]
        data_temp['DNAFISH_mean_int_merge'] = [DNAFISH_props[i].intensity_mean for i in range(len(DNAFISH_props))]
        data_temp['bg_merge'] = [bg_merge] * len(label_props)

        data = pd.concat([data, data_temp], axis=0)

data['DNAFISH_total_int_merge'] = [(data['DNAFISH_mean_int_merge'].tolist()[i]-data['bg_merge'].tolist()[i]) *
                                   data['nuclear_area_merge'].tolist()[i] if data['DNAFISH_mean_int_merge'].tolist()[i]-data['bg_merge'].tolist()[i] >0
                                   else 0 for i in range(len(data))]

data.to_csv('%s/%s/16_%s_copy_number.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 17-1
print("Running 17-1...")
local_size = 200
data = pd.DataFrame(columns=['sample', 'fov', 'DNAFISH_thresh'])
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_DNAFISH = img_stack[fov, :, 1, :, :]
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_hoechst_seg_z = img_hoechst[seg_z]
        img_DNAFISH_seg_z = img_DNAFISH[seg_z]

        img_seg_z = skio.imread("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_ecDNA = np.zeros_like(img_DNAFISH_seg_z)
        thresh = threshold_otsu(img_DNAFISH_seg_z)
        img_seg_ecDNA[img_DNAFISH_seg_z > thresh] = 1
        data.loc[len(data.index)] = [s, fov, thresh]

        if not os.path.exists("%s/%s/17_seg_DNAFISH_tif/" % (output_dir, sample)):
            os.makedirs("%s/%s/17_seg_DNAFISH_tif/" % (output_dir, sample))
        tif.imwrite("%s/%s/17_seg_DNAFISH_tif/%s_%s_seg_DNAFISH.tif" % (output_dir, sample, s, fov), img_seg_ecDNA)

        viewer = napari.Viewer()
        viewer.add_image(img_DNAFISH_seg_z, blending='additive', colormap='green', contrast_limits=[0, 20000])
        viewer.add_image(img_seg_ecDNA, blending='additive', contrast_limits=[0, 3])
        if not os.path.exists("%s/%s/17_seg_DNAFISH_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/17_seg_DNAFISH_color/" % (output_dir, sample))
        plt.imsave("%s%s/17_seg_DNAFISH_color/%s_%s_seg_DNAFISH_color.tiff" % (output_dir, sample, s, fov),
                   dis.blending(viewer))
        viewer.close()

data.to_csv('%s/%s/17_%s_DNAFISH_thresh.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 17-2
print("Running 17-2...")
data = pd.DataFrame(columns=['sample', 'fov', 'label_nuclear', 'area_nuclear',
                             'n_ecDNA', 'area_ind_ecDNA', 'total_area_ecDNA',
                             'area_ratio_ind_ecDNA', 'total_area_ratio_ecDNA',
                             'percentage_area_curve_ecDNA', 'percentage_area_n_half',
                             'cum_area_ind_ecDNA', 'cum_area_n_half', 'per_AUC'])

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
                per_AUC = np.sum(cum_percentage_area_ind_ecDNA[:10])/10
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

### 18
print("Running 18...")
data = pd.DataFrame(columns=['sample', 'fov', 'label_nuclear',
                             'DNAFISH_seg_label', 'int_r_to_edge', 'int_relative_r'])


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov + 1, total_fov))
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
            local_nuclear_props = regionprops(label(local_nuclear_seg), local_nuclear)

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

            data.loc[len(data.index)] = [s, fov, label_nuclear, DNAFISH_seg_label, int_r_to_edge, int_relative_r]

data.to_csv('%s/%s/18_%s_radial.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")
