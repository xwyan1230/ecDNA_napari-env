import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_DNAFISH_ctrl/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230409_DM_plate_DNAFISH_control'
sample = 'B2'
batch = 0
file_name = '%s_%s_RAW' % (exp, sample)
img_hoechst = skio.imread("%s/%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s/%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")

print(img_hoechst.shape)

for i in range(img_hoechst.shape[0]):
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst[i, :, :], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH[i, :, :], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    napari.run()

