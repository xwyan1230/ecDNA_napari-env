import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_H4K16Ac/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'HSR'
total_fov = 49
start_fov = 0
for f in range(total_fov):
    fov = start_fov + f
    print(fov)
    if fov < 10:
        file_name = '20230721_H4K16Ac_%s_H4K16Ac_s0%s' % (sample, fov)
    else:
        file_name = '20230721_H4K16Ac_%s_H4K16Ac_s%s' % (sample, fov)
    img_hoechst = skio.imread("%sraw/%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%sraw/%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_IF = skio.imread("%sraw/%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
    # img_hoechst_temp = np.concatenate([np.zeros(shape=[3, 3144], dtype=np.uint16), img_hoechst_temp], axis=0)[:3144, :3144]
    # img_DNAFISH_temp = np.concatenate([np.zeros(shape=[6, 3144], dtype=np.uint16), img_DNAFISH_temp], axis=0)[:3144, :3144]

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_IF, blending='additive', colormap='red', contrast_limits=[0, img_IF.max()])
    napari.run()