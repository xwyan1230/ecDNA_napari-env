import napari
import skimage.io as skio
import matplotlib.pyplot as plt
import os
import numpy as np

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
sample = 'D3'
raw_folder = '01_raw'
save_folder = 'v5_raw-image'
save_path = '%s%s/' % (master_folder, save_folder)
if not os.path.exists(save_path):
    os.makedirs(save_path)

fov = 2

im_z_stack_nuclear = skio.imread("%s%s/%s/%s/R%s_RAW_ch00.tif" %
                                 (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")
im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s/R%s_RAW_ch01.tif" %
                                 (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")
im_z_stack_IF = skio.imread("%s%s/%s/%s/R%s_RAW_ch02.tif" %
                            (master_folder, sample[0], sample[1:], raw_folder, fov), plugin="tifffile")

viewer = napari.view_image(im_z_stack_nuclear, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(im_z_stack_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 8000])
viewer.add_image(im_z_stack_IF, blending='additive', colormap='magenta', contrast_limits=[0, 20000])
napari.run()