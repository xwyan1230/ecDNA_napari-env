import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import math
import nd2
import shared.segmentation as seg
import shared.objects as obj
from skimage.measure import label, regionprops
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'D4'
samples = ['D4']
total_fovs = [16]

### 09
convex_conversion_threshold = 0.8
circ_threshold = 0.8
local_size = 200
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = 3000
max_size_nuclear = 9000

data = pd.DataFrame(columns=['sample', 'fov', 'seg_z'])

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)

        img_nuclear_seg = seg.nuclear_seg_nikon(img_hoechst_merge, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)
        img_nuclear_seg_convex = obj.label_remove_low_circ(seg.obj_to_convex_filter(img_nuclear_seg, threshold=convex_conversion_threshold), thresh=circ_threshold)

        if not os.path.exists("%s/%s/09_seg_tif/" % (output_dir, sample)):
            os.makedirs("%s/%s/09_seg_tif/" % (output_dir, sample))
        if not os.path.exists("%s/%s/09_seg_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/09_seg_color/" % (output_dir, sample))
        tif.imwrite("%s/%s/09_seg_tif/%s_%s_seg.tif" % (output_dir, sample, s, fov), img_nuclear_seg_convex)

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
        # viewer.add_image(img_nuclear_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
        viewer.add_image(img_nuclear_seg_convex, blending='additive', contrast_limits=[0, 1])
        plt.imsave("%s%s/09_seg_color/%s_%s_seg_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

        ### 10

        total_hoechst_lst = []
        for z in range(img_hoechst.shape[0]):
            total_hoechst = np.sum(img_hoechst[z])
            total_hoechst_lst.append(total_hoechst)
        seg_z = total_hoechst_lst.index(max(total_hoechst_lst))
        print(seg_z)
        data.loc[len(data.index)] = [s, fov, seg_z]
        img_hoechst_seg = img_hoechst[seg_z]

        img_nuclear_seg = seg.nuclear_seg_nikon(img_hoechst_seg, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                                max_size=max_size_nuclear)
        img_nuclear_seg_convex = obj.label_remove_low_circ(seg.obj_to_convex_filter(img_nuclear_seg, threshold=convex_conversion_threshold), thresh=circ_threshold)

        if not os.path.exists("%s/%s/10_seg_z_tif/" % (output_dir, sample)):
            os.makedirs("%s/%s/10_seg_z_tif/" % (output_dir, sample))
        if not os.path.exists("%s/%s/10_seg_z_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/10_seg_z_color/" % (output_dir, sample))
        tif.imwrite("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, s, fov), img_nuclear_seg_convex)

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_seg, blending='additive', colormap='blue', contrast_limits=[0, 12000])
        # viewer.add_image(img_nuclear_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
        viewer.add_image(img_nuclear_seg_convex, blending='additive', contrast_limits=[0, 1])
        plt.imsave("%s%s/10_seg_z_color/%s_%s_seg_z_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

data.to_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 11
data = pd.DataFrame(columns=['sample', 'fov', 'label', 'seg', 'centroid_0', 'centroid_1', 'mean_int'])

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_seg_total = skio.imread("%s/%s/09_seg_tif/%s_%s_seg.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_z = skio.imread("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, s, fov), plugin="tifffile")

        total_props = regionprops(label(img_seg_total), img_seg_total)
        for i in range(len(total_props)):
            data.loc[len(data.index)] = [s, fov, total_props[i].label, 'total', total_props[i].centroid[0], total_props[i].centroid[1], total_props[i].intensity_mean]
        z_props = regionprops(label(img_seg_z), img_seg_z)
        for i in range(len(z_props)):
            data.loc[len(data.index)] = [s, fov, z_props[i].label, 'z', z_props[i].centroid[0], z_props[i].centroid[1], z_props[i].intensity_mean]

data.to_csv('%s/%s/11_%s_seg_centroids.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 12
data_common = pd.DataFrame(columns=['sample', 'fov', 'label', 'seg', 'centroid_0', 'centroid_1', 'mean_int', 'label_new'])

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        data_total = data[(data['sample'] == s) & (data['fov'] == fov) & (data['seg'] == 'total')].copy().reset_index(drop=True)
        data_z = data[(data['sample'] == s) & (data['fov'] == fov) & (data['seg'] == 'z')].copy().reset_index(drop=True)
        count = 1
        for i in range(len(data_total)):
            print("%s/%s" % (i+1, len(data_total)))
            for j in range(len(data_z)):
                if (data_z['centroid_0'][j]>(data_total['centroid_0'][i]-5)) & (data_z['centroid_0'][j]<(data_total['centroid_0'][i]+5)) & (data_z['centroid_1'][j]>(data_total['centroid_1'][i]-5)) & (data_z['centroid_1'][j]<(data_total['centroid_1'][i]+5)):
                    data_common.loc[len(data_common.index)] = [s, fov, data_total['label'][i], 'total', data_total['centroid_0'][i], data_total['centroid_1'][i], data_total['mean_int'][i], count]
                    data_common.loc[len(data_common.index)] = [s, fov, data_z['label'][j], 'z', data_z['centroid_0'][j], data_z['centroid_1'][j], data_z['mean_int'][j], count]
                    count = count + 1
                    break

data_common.to_csv('%s/%s/12_%s_seg_common_centroids.txt' % (output_dir, sample, sample), index=False, sep='\t')

### 13
data = pd.read_csv("%s/%s/12_%s_seg_common_centroids.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')
pd_seg_z = pd.read_csv('%s/%s/10_%s_seg_z.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')

for k in range(len(samples)):
    s = samples[k]
    total_fov = total_fovs[k]
    img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, s))
    pd_seg_z_s = pd_seg_z[pd_seg_z['sample'] == s].copy().reset_index(drop=True)

    for fov in range(total_fov):
        print("%s/%s" % (fov+1, total_fov))
        img_hoechst = img_stack[fov, :, 0, :, :]
        img_hoechst_merge = img_hoechst.max(axis=0)
        seg_z = pd_seg_z_s['seg_z'][fov]
        img_hoechst_seg_z = img_hoechst[seg_z]
        img_seg_total = skio.imread("%s/%s/09_seg_tif/%s_%s_seg.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_z = skio.imread("%s/%s/10_seg_z_tif/%s_%s_seg_z.tif" % (output_dir, sample, s, fov), plugin="tifffile")
        img_seg_total_new = np.zeros_like(img_seg_total)
        img_seg_z_new = np.zeros_like(img_seg_z)
        data_total = data[(data['sample'] == s) & (data['fov'] == fov) & (data['seg'] == 'total')].copy().reset_index(drop=True)
        data_z = data[(data['sample'] == s) & (data['fov'] == fov) & (data['seg'] == 'z')].copy().reset_index(drop=True)

        total_props = regionprops(label(img_seg_total), img_seg_total)
        for i in range(len(total_props)):
            if total_props[i].intensity_mean in data_total['mean_int'].tolist():
                label_index = data_total[data_total['mean_int'] == total_props[i].intensity_mean].index
                print('total: %s-%s' % (total_props[i].intensity_mean, data_total['label_new'][label_index]))
                img_seg_total_new[img_seg_total == total_props[i].intensity_mean] = data_total['label_new'][label_index]

        z_props = regionprops(label(img_seg_z), img_seg_z)
        for i in range(len(z_props)):
            if z_props[i].intensity_mean in data_z['mean_int'].tolist():
                label_index = data_z[data_z['mean_int'] == z_props[i].intensity_mean].index
                print('z: %s-%s' % (z_props[i].intensity_mean, data_z['label_new'][label_index]))
                img_seg_z_new[img_seg_z == z_props[i].intensity_mean] = data_z['label_new'][label_index]

        if not os.path.exists("%s/%s/13_seg_new_tif/" % (output_dir, sample)):
            os.makedirs("%s/%s/13_seg_new_tif/" % (output_dir, sample))
        tif.imwrite("%s/%s/13_seg_new_tif/%s_%s_seg_new.tif" % (output_dir, sample, s, fov), img_seg_total_new)
        tif.imwrite("%s/%s/13_seg_new_tif/%s_%s_seg_z_new.tif" % (output_dir, sample, s, fov), img_seg_z_new)

        if not os.path.exists("%s/%s/13_seg_new_color/" % (output_dir, sample)):
            os.makedirs("%s/%s/13_seg_new_color/" % (output_dir, sample))

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_merge, blending='additive', colormap='blue', contrast_limits=[0, 12000])
        viewer.add_image(img_seg_total_new, blending='additive', contrast_limits=[0, 1])
        plt.imsave("%s%s/13_seg_new_color/%s_%s_seg_new_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_seg_z, blending='additive', colormap='blue', contrast_limits=[0, 12000])
        viewer.add_image(img_seg_z_new, blending='additive', contrast_limits=[0, 1])
        plt.imsave("%s%s/13_seg_new_color/%s_%s_seg_z_new_color.tiff" % (output_dir, sample, s, fov), dis.blending(viewer))
        viewer.close()

print("DONE!")