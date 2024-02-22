import skimage.io as skio
import napari
import imutils
import shared.image as ima
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import shared.display as dis
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230711_analysis_H3K4me3_NPC/"
data_dir = "%sdata/raw/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM'
# sample = 'NPC'
sample = 'TPR'
# prefix = '20230616_H3K4me3_DM_H3K4me3'
# prefix = '20230608_DM_IF-NPC_acidFISH-MYC_ColoDM_IF-NPC_acidFISH-MYC_100x-2048'
prefix = '20230616_DM_TPR-Nup98_DM_TPR'
total_fov = 2
start_fov = 1
fov = 1
z = 13
center = [1543.9267693, 904.96859135]
local_range = 200
position = [int(center[0]-local_range), int(center[0]+local_range), int(center[1]-local_range), int(center[1]+local_range)]

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.3, 0.05679, 0.05679])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]

lst = [x.split('%s_' % prefix)[1].split('_ch')[0] for x in os.listdir('%s%s' % (data_dir, sample)) if x[-4:] == '.tif']

lst_temp = [int(x.split('_z')[1]) for x in lst if x.split('_z')[0] == '%s' % fov]
total_z = max(lst_temp) + 1
img_3d_nuclear = skio.imread("%s%s/%s_%s_z%s_ch00.tif" % (data_dir, sample, prefix, fov, z), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]
img_3d_DNAFISH = skio.imread("%s%s/%s_%s_z%s_ch01.tif" % (data_dir, sample, prefix, fov, z), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]
img_3d_laminB = skio.imread("%s%s/%s_%s_z%s_ch02.tif" % (data_dir, sample, prefix, fov, z), plugin="tifffile")[np.newaxis, position[0]:position[1], position[2]:position[3]]

if not os.path.exists('%s/color_img/%s/fov%s_%s_%s_tif/' % (output_dir, sample, fov, center[0], center[1])):
    os.makedirs('%s/color_img/%s/fov%s_%s_%s_tif/' % (output_dir, sample, fov, center[0], center[1]))
tif.imwrite('%s/color_img/%s/fov%s_%s_%s_tif/z%s_nuclear.tiff' % (output_dir, sample, fov, center[0], center[1], z), img_3d_nuclear)
tif.imwrite('%s/color_img/%s/fov%s_%s_%s_tif/z%s_DNAFISH.tiff' % (output_dir, sample, fov, center[0], center[1], z), img_3d_DNAFISH)
tif.imwrite('%s/color_img/%s/fov%s_%s_%s_tif/z%s_%s.tiff' % (output_dir, sample, fov, center[0], center[1], z, sample), img_3d_laminB)

print("DONE!")

