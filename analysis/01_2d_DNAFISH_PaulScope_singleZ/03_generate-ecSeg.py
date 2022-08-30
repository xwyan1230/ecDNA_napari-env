from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion, binary_erosion, binary_dilation, extrema, disk
import pandas as pd
import numpy as np
import shared.image as ima
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg
import tifffile as tif
import shared.math as mat
from skimage import segmentation
import shared.display as dis
from skimage.filters import threshold_otsu, threshold_local, sobel
import shared.objects as obj
import matplotlib.pyplot as plt
import os
import napari

# INPUT PARAMETERS
# file info
# master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220816_Natasha_ColoDM_reimage/"
master_folder = '/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220825_POLR3D/20220825_POLR3Dtest/'
sample = 'Control_2'
master_path = '%s%s/' % (master_folder, sample)
threshold_max = 25000
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100
rmax = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s_z.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

fovs = list(set(data_z['FOV'].tolist()))
total_fov = len(fovs)

for f in range(total_fov):
    print("Analyzing %s, analyzing fov %s/%s" % (sample, f+1, total_fov))
    fov = fovs[f]
    data_z_fov = data_z[data_z['FOV'] == fov].copy().reset_index()

    # load images
    file_prefix = "%s_Position %s_RAW" % (sample, fov)
    im_z_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s_ch01.tif" % (master_path, file_prefix), plugin="tifffile")
    im_z_stack_seg_convex = skio.imread("%s%s_seg.tif" % (master_path, file_prefix), plugin="tifffile")

    im_z_stack_nuclear_seg_z = np.zeros(shape=im_z_stack_nuclear.shape, dtype=np.uint16)
    im_z_stack_DNAFISH_seg_z = np.zeros(shape=im_z_stack_nuclear.shape, dtype=np.uint16)
    im_z_stack_nuclear_max = np.zeros(shape=im_z_stack_nuclear.shape, dtype=np.uint16)
    im_z_stack_DNAFISH_max = np.zeros(shape=im_z_stack_nuclear.shape, dtype=np.uint16)

    for i in range(len(data_z_fov)):
        print("Analyzing %s, analyzing fov %s, nucleus %s/%s" % (sample, f+1, i+1, len(data_z_fov)))
        z_current = data_z_fov['z'][i]
        label_nuclear = data_z_fov['label_nuclear'][i]
        original_centroid_nuclear = data_z_fov['centroid_nuclear'][i]

        # get images for given z
        img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
        img_nuclear = im_z_stack_nuclear[z_current]
        img_DNAFISH = im_z_stack_DNAFISH[z_current]

        # get local images
        position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_nuclear_centroid = local_nuclear_props[0].centroid

        # ecDNA segmentation
        local_DNAFISH_singlet = local_DNAFISH.copy()
        local_DNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
        otsu_threshold_val_local_DNAFISH = threshold_otsu(local_DNAFISH_singlet)
        threshold = otsu_threshold_val_local_DNAFISH + (threshold_max - otsu_threshold_val_local_DNAFISH) / 10

        k_dots = 5000
        vector = []
        vector_cum_weight = []
        weight = 0

        DNAFISH_seg, _ = seg.find_organelle(local_DNAFISH, 'na', extreme_val=500, bg_val=data_z['limit'][i] * 0.8,
                                            min_size=0, max_size=500)
        DNAFISH_seg[local_nuclear_seg_convex == 0] = 0
        DNAFISH_seg = binary_dilation(DNAFISH_seg)

        bg_seg_temp = local_nuclear_seg_convex.copy()
        bg_seg_temp[DNAFISH_seg == 1] = 0
        local_DNAFISH_bg_temp = local_DNAFISH.copy()
        local_DNAFISH_bg_temp[local_nuclear_seg_convex == 0] = 0
        local_DNAFISH_bg_temp[DNAFISH_seg == 1] = 0
        bg = np.sum(local_DNAFISH_bg_temp) / np.sum(bg_seg_temp)

        int_thresh = bg * 1.2
        for m in range(len(local_nuclear_seg_convex)):
            for n in range(len(local_nuclear_seg_convex[0])):
                if local_nuclear_seg_convex[m][n] == 1:
                    vector.append([m, n])
                    if local_DNAFISH[m][n] > int_thresh:
                        weight = weight + local_DNAFISH[m][n] - int_thresh
                    vector_cum_weight.append(weight)
        if weight != 0:
            random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=int(k_dots * weight / 50000))
            img_dot = np.zeros_like(local_nuclear_seg_convex)
            for m in random_dot:
                img_dot[m[0]][m[1]] = img_dot[m[0]][m[1]] + 1
            img_dot_remove_bg = img_dot.copy()
            img_dot_remove_bg[img_dot < 35] = 0
            img_dot_remove_bg_thresh = img_dot_remove_bg > threshold_otsu(img_dot_remove_bg)
            img_dot_remove_bg_thresh = obj.remove_small(label(img_dot_remove_bg_thresh), 5)
            img_dot_remove_bg_thresh = erosion(img_dot_remove_bg_thresh)
            img_dot_remove_bg_thresh = erosion(img_dot_remove_bg_thresh)
            img_dot_remove_bg_thresh = dilation(img_dot_remove_bg_thresh)
            # img_dot_remove_bg_thresh = dilation(img_dot_remove_bg_thresh)

            ecDNA_seg = img_dot_remove_bg_thresh.copy()
            ecDNA_seg[DNAFISH_seg == 1] = 1
            ecDNA_seg_filter = np.zeros_like(ecDNA_seg)
            ecDNA_seg_label = label(ecDNA_seg)
            ecDNA_seg_props = regionprops(ecDNA_seg_label, local_DNAFISH)
            for j in range(len(ecDNA_seg_props)):
                if ecDNA_seg_props[j].intensity_mean > threshold:
                    ecDNA_seg_filter[ecDNA_seg_label == ecDNA_seg_props[j].label] = 1

            im_z_stack_DNAFISH_seg_z[z_current] = \
                ima.image_paste_to(im_z_stack_DNAFISH_seg_z[z_current], ecDNA_seg_filter,
                                   [int(original_centroid_nuclear[0] - local_nuclear_centroid[0]),
                                    int(original_centroid_nuclear[1] - local_nuclear_centroid[1])])
            im_z_stack_nuclear_seg_z[z_current] = \
                ima.image_paste_to(im_z_stack_nuclear_seg_z[z_current], local_nuclear_seg_convex,
                                   [int(original_centroid_nuclear[0] - local_nuclear_centroid[0]),
                                    int(original_centroid_nuclear[1] - local_nuclear_centroid[1])])
            im_z_stack_DNAFISH_max[z_current] = im_z_stack_DNAFISH[z_current]
            im_z_stack_nuclear_max[z_current] = im_z_stack_nuclear[z_current]

    tif.imwrite("%s%s_ecSeg_z.tif" % (master_path, file_prefix), im_z_stack_DNAFISH_seg_z)
    tif.imwrite("%s%s_seg_z.tif" % (master_path, file_prefix), im_z_stack_nuclear_seg_z)

    img_nuclear_seg_z_max = im_z_stack_nuclear_seg_z.max(axis=0)
    img_DNAFISH_seg_z_max = im_z_stack_DNAFISH_seg_z.max(axis=0)
    img_nuclear_max = im_z_stack_nuclear_max.max(axis=0)
    img_DNAFISH_max = im_z_stack_DNAFISH_max.max(axis=0)

    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH_max, blending='additive', colormap='green', contrast_limits=[0, threshold_max])
    viewer.add_image(img_nuclear_max, blending='additive', colormap='blue')
    plt.imsave('%s%s_napari-img.tiff' % (master_path, file_prefix), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH_seg_z_max, blending='additive', colormap='green')
    viewer.add_image(img_nuclear_seg_z_max, blending='additive', colormap='blue')
    plt.imsave('%s%s_napari-seg.tiff' % (master_path, file_prefix), dis.blending(viewer))
    viewer.close()

    """save_folder = '%snapari_z/' % master_path
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    for k in range(im_z_stack_nuclear_seg_z.shape[0]):
        if im_z_stack_nuclear_seg_z[k].max() != 0:
            viewer = napari.Viewer()
            viewer.add_image(im_z_stack_DNAFISH[k], blending='additive', colormap='green', contrast_limits=[0, threshold_max])
            viewer.add_image(im_z_stack_nuclear[k], blending='additive', colormap='blue')
            plt.imsave('%s%s_napari-img_z%s.tiff' % (save_folder, file_prefix, k), dis.blending(viewer))
            viewer.close()

            viewer = napari.Viewer()
            viewer.add_image(im_z_stack_DNAFISH_seg_z[k], blending='additive', colormap='green')
            viewer.add_image(im_z_stack_nuclear_seg_z[k], blending='additive', colormap='blue')
            plt.imsave('%s%s_napari-seg_z%s.tiff' % (save_folder, file_prefix, k), dis.blending(viewer))
            viewer.close()"""

print("DONE!")
