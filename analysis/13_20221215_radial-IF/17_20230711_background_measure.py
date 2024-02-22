import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H3K4me3'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s%s_ch02.tif" % (data_dir, file_name), plugin="tifffile")
img_IF = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")

print(img_hoechst.shape)

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        img_hoechst_temp = img_hoechst[i, :, :]
        img_DNAFISH_temp = img_DNAFISH[i, :, :]
        img_IF_temp = img_IF[i, :, :]
        img_hoechst_temp = np.concatenate([np.zeros(shape=[3, 3144], dtype=np.uint16), img_hoechst_temp], axis=0)[:3144, :3144]
        img_DNAFISH_temp = np.concatenate([np.zeros(shape=[6, 3144], dtype=np.uint16), img_DNAFISH_temp], axis=0)[:3144, :3144]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst_temp[:3144, :3144], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH_temp[:3144, :3144], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        viewer.add_image(img_IF_temp[:3144, :3144], blending='additive', colormap='red', contrast_limits=[0, img_IF.max()])
        napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_IF, blending='additive', colormap='green', contrast_limits=[0, img_IF.max()])
    napari.run()
