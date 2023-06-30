import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/data/DM_DNAFISH/"
data_dir = "%s" % master_folder
output_dir = "%s" % master_folder

sample = 'DM_3_49pos'
file_name = '%s_RAW' % sample
img_hoechst = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")

print(img_hoechst.shape)

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, :, :], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH[i, :, :], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_DNAFISH.max()])
    napari.run()
