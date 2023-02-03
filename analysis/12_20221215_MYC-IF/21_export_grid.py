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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-old'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_MYC_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")
hue_order = [10000, 15000, 20000, 25000, 30000, 35000, 40000, 70000]

n_nuclear_convex_dilation = 8
local_size = 200
rmax = 100
size = 20*local_size + 22

group_num = 8
cell_count = [0] * group_num
img_grid_hoechst = np.zeros(shape=(group_num, size, size), dtype=np.uint16)
img_grid_DNAFISH = np.zeros(shape=(group_num, size, size), dtype=np.uint16)
img_grid_MYCIF = np.zeros(shape=(group_num, size, size), dtype=np.uint16)


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


for fov in range(img_MYC_stack.shape[0]):
    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_MYC = img_MYC_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s%s/seg_tif/%s_%s_ecseg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    # measure
    # get local images
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):

        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_mCherry = img_MYC.copy()
        local_mCherry = local_mCherry[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_nuclear_seg_mCherry = ima.img_local_seg(img_seg, position, nuclear_props[i].label)
        local_nuclear[local_nuclear_seg == 0] = 0
        local_DNAFISH[local_nuclear_seg == 0] = 0
        local_mCherry[local_nuclear_seg == 0] = 0

        # basic measurements
        local_IF_props = regionprops(label(local_nuclear_seg_mCherry), local_mCherry)
        mean_int_IF = local_IF_props[0].intensity_mean
        local_nuclear_props = regionprops(label(local_nuclear_seg))
        local_centroid = local_nuclear_props[0].centroid
        print(mean_int_IF)
        pos = dat.find_pos(mean_int_IF, hue_order)
        print(pos)
        paste_to_centroid = [grid_location(cell_count[pos])[0], grid_location(cell_count[pos])[1]]
        print(paste_to_centroid)
        cell_count[pos] += 1
        img_grid_hoechst[pos] = ima.image_paste_to(img_grid_hoechst[pos], local_nuclear, [int(paste_to_centroid[0] - local_centroid[0]), int(paste_to_centroid[1] - local_centroid[1])])
        img_grid_DNAFISH[pos] = ima.image_paste_to(img_grid_DNAFISH[pos], local_DNAFISH, [int(paste_to_centroid[0] - local_centroid[0]), int(paste_to_centroid[1] - local_centroid[1])])
        img_grid_MYCIF[pos] = ima.image_paste_to(img_grid_MYCIF[pos], local_mCherry, [int(paste_to_centroid[0] - local_centroid[0]), int(paste_to_centroid[1] - local_centroid[1])])

folder = 'grid_seg'
if not os.path.exists("%s%s/%s/" % (output_dir, sample, folder)):
    os.makedirs("%s%s/%s/" % (output_dir, sample, folder))
tif.imwrite("%s%s/%s/%s_hoechst.tif" % (output_dir, sample, folder, sample), img_grid_hoechst)
tif.imwrite("%s%s/%s/%s_DNAFISH.tif" % (output_dir, sample, folder, sample), img_grid_DNAFISH)
tif.imwrite("%s%s/%s/%s_MYCIF.tif" % (output_dir, sample, folder, sample), img_grid_MYCIF)

for i in range(group_num):
    viewer = napari.Viewer()
    viewer.add_image(img_grid_hoechst[i], blending='additive', colormap='blue', contrast_limits=[0, img_grid_hoechst[i].max()])
    viewer.add_image(img_grid_DNAFISH[i], blending='additive', colormap='green', contrast_limits=[0, img_grid_DNAFISH[i].max()])
    viewer.add_image(img_grid_MYCIF[i], blending='additive', colormap='magenta', contrast_limits=[0, img_grid_MYCIF[i].max()])
    plt.imsave('%s%s/%s/%s_group%s.tiff' % (output_dir, sample, folder, sample, hue_order[i]), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_grid_hoechst[i], blending='additive', colormap='blue',
                     contrast_limits=[0, img_grid_hoechst[i].max()])
    viewer.add_image(img_grid_DNAFISH[i], blending='additive', colormap='green',
                     contrast_limits=[0, img_grid_DNAFISH[i].max()])
    plt.imsave('%s%s/%s/%s_group%s_DNAFISH.tiff' % (output_dir, sample, folder, sample, hue_order[i]), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_grid_MYCIF[i], blending='additive', colormap='magenta',
                     contrast_limits=[0, img_grid_MYCIF[i].max()])
    plt.imsave('%s%s/%s/%s_group%s_MYC.tiff' % (output_dir, sample, folder, sample, hue_order[i]), dis.blending(viewer))
    viewer.close()