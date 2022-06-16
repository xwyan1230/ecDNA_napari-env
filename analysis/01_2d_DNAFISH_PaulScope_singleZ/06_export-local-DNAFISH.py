from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import numpy as np
import shared.image as img
import skimage.io as skio
import shared.dataframe as dat
import random
import shared.segmentation as seg
import tifffile as tif
import matplotlib.pyplot as plt
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = '72hrJQ1'
raw_folder = '01_raw'
seg_folder = '02_seg'
save_folder = '04_DNAFISH_scale'
total_fov = 18
total_z = 23
# cell info
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, Paul scope 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
z_size = 500  # nm (Paul scope)
local_size = 100

# LOAD Z FILE
data_z = pd.read_csv('%s%s/%s_z.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

# IMAGING ANALYSIS
for i in range(len(data_z)):
    print("Analyzing nucleus %s/%s" % (i+1, len(data_z)))
    fov = data_z['FOV'][i]
    z_current = data_z['z'][i]
    label_nuclear = data_z['label_nuclear'][i]
    original_centroid_nuclear = data_z['centroid_nuclear'][i]
    # load images
    im_z_stack_nuclear = skio.imread("%s%s/%s/%s_RAW_ch00_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s_RAW_ch01_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    im_z_stack_seg_convex = skio.imread("%s%s/%s/%s_seg_fov%s.tif" % (master_folder, sample, seg_folder, sample, fov),
                                        plugin="tifffile")
    # get images for given z
    img_nuclear_seg_convex = im_z_stack_seg_convex[z_current]
    img_nuclear = im_z_stack_nuclear[z_current]
    img_DNAFISH = im_z_stack_DNAFISH[z_current]
    # get local images
    position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
    local_nuclear = img_nuclear.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]

    # ecDNA segmentation
    int_thresh = data_z['limit'][i]*0.9
    k_dots = 5000
    vector = []
    vector_cum_weight = []
    weight = 0
    for m in range(len(local_nuclear_seg_convex)):
        for n in range(len(local_nuclear_seg_convex[0])):
            if local_nuclear_seg_convex[m][n] == 1:
                vector.append([m, n])
                if local_DNAFISH[m][n] > int_thresh:
                    weight = weight + local_DNAFISH[m][n] - int_thresh
                vector_cum_weight.append(weight)
    if weight != 0:
        random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=k_dots)
        img_dot = np.zeros_like(local_nuclear_seg_convex)
        for m in random_dot:
            img_dot[m[0]][m[1]] = img_dot[m[0]][m[1]] + 1
        img_dot_remove_bg = dilation(img_dot)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = dilation(img_dot_remove_bg)
        img_dot_seg = img_dot_remove_bg.copy()
        img_dot_seg[img_dot_remove_bg > 0] = 1
        DNAFISH_seg, _ = seg.find_organelle(local_DNAFISH, 'na', extreme_val=500, bg_val=data_z['limit'][i]*0.8, min_size=0,
                                            max_size=10000)
        DNAFISH_seg[local_nuclear_seg_convex == 0] = 0
        ecDNA_seg = img_dot_seg.copy()
        ecDNA_seg[DNAFISH_seg == 1] = 1
        plt.imsave('%s%s/%s/%s_DNAFISH_fov%s_z%s_i%s.tiff' % (master_folder, sample, save_folder, sample, fov, z_current,
                                                              label_nuclear), local_DNAFISH, vmin=0, vmax=6000)
        plt.imsave('%s%s/%s/%s_ecSeg_fov%s_z%s_i%s.tiff' % (master_folder, sample, save_folder, sample, fov, z_current,
                                                            label_nuclear), ecDNA_seg, vmin=0, vmax=1)

print("DONE!")
