import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import tifffile as tif
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


def grid_location(num: int):
    row = int(num / 10)
    column = num % 10
    x = row * (2 * local_size + 2) + 2 + local_size
    y = column * (2 * local_size + 2) + 2 + local_size
    return x, y


n_nuclear_convex_dilation = 0
local_size = 150
size = 20*local_size + 22

sample = 'PC9'
n = 0
# fov_lst = [fov_to_str(i) for i in np.arange(2, 37, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 7, 1)]  # Colo320DM
# fov_lst = [fov_to_str(i) for i in np.arange(1, 31, 1)]  # Colo320HSR
# fov_lst = ['MycI2_'+fov_to_str(i) for i in np.arange(1, 17, 1)]  # HCT116
# fov_lst = list(np.arange(1, 9, 1)) + list(np.arange(10, 13, 1)) + [fov_to_str(i) for i in np.arange(13, 16, 1)] + [fov_to_str(i) for i in np.arange(17, 26, 1)] + [fov_to_str(i) for i in np.arange(27, 36, 1)] # PC3
fov_lst = [fov_to_str(i) for i in np.arange(1, 2, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 8, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(9, 19, 1)]  # PC9

img_count = 0
cell_count = 0
img_grid_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_grid_RNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

for fov in range(len(fov_lst)):

    if cell_count == 100:
        if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
            os.makedirs("%s%s/grid/" % (output_dir, sample))
        tif.imwrite("%s%s/grid/%s_hoechst_%s.tif" % (output_dir, sample, sample, img_count), img_grid_hoechst)
        tif.imwrite("%s%s/grid/%s_RNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_grid_RNAFISH)

        viewer = napari.Viewer()
        viewer.add_image(img_grid_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_grid_hoechst.max()])
        viewer.add_image(img_grid_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_grid_RNAFISH.max()])
        plt.imsave('%s%s/grid/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
        viewer.close()

        cell_count = 0
        img_count += 1
        img_grid_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_grid_RNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

    print(fov)

    file_name = '%s_%s_Lng_SVCC_Processed001_RAW' % (sample, fov_lst[fov])
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_RNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir1, sample, file_name), plugin="tifffile")
    img_seg = skio.imread("%s%s/seg_tif_manual/%s_%s_seg.tif" % (data_dir2, sample, sample, fov_lst[fov]), plugin="tifffile")

    nuclear_props = regionprops(img_seg)

    for i in range(len(nuclear_props)):
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_seg, position, nuclear_props[i].label)
        local_nuclear = img_hoechst.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_RNAFISH = img_RNAFISH.copy()
        local_RNAFISH = local_RNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid
        local_seg_filter = local_nuclear_seg_convex
        local_nuclear[local_seg_filter == 0] = 0
        local_RNAFISH[local_seg_filter == 0] = 0

        if cell_count < 100:
            paste_to_centroid = [grid_location(cell_count)[0], grid_location(cell_count)[1]]
            print(paste_to_centroid)
            cell_count += 1
            img_grid_hoechst = ima.image_paste_to(img_grid_hoechst, local_nuclear, [int(paste_to_centroid[0] - local_centroid[0]), int(paste_to_centroid[1] - local_centroid[1])])
            img_grid_RNAFISH = ima.image_paste_to(img_grid_RNAFISH, local_RNAFISH, [int(paste_to_centroid[0] - local_centroid[0]), int(paste_to_centroid[1] - local_centroid[1])])

if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
    os.makedirs("%s%s/grid/" % (output_dir, sample))
tif.imwrite("%s%s/grid/%s_hoechst_%s.tif" % (output_dir, sample, sample, img_count), img_grid_hoechst)
tif.imwrite("%s%s/grid/%s_RNAFISH_%s.tif" % (output_dir, sample, sample, img_count), img_grid_RNAFISH)

viewer = napari.Viewer()
viewer.add_image(img_grid_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_grid_hoechst.max()])
viewer.add_image(img_grid_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_grid_RNAFISH.max()])
plt.imsave('%s%s/grid/%s_%s.tiff' % (output_dir, sample, sample, img_count), dis.blending(viewer))
viewer.close()