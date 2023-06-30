import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230503_analysis_Natasha_TPL/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'TPL_mCh-DMSO'
batch = 1
start_fov = 2
end_fov = 51
total_fov = end_fov-start_fov+1

for f in range(total_fov):
    fov = start_fov+f
    file_name = 'Position %s' % fov
    img_hoechst = skio.imread("%s%s/%s/%s/%s_ch00.tif" % (data_dir, sample, batch, file_name, file_name), plugin="tifffile")
    img_DNAFISH = skio.imread("%s%s/%s/%s/%s_ch02.tif" % (data_dir, sample, batch, file_name, file_name), plugin="tifffile")
    img_mCherry = skio.imread("%s%s/%s/%s/%s_ch01.tif" % (data_dir, sample, batch, file_name, file_name), plugin="tifffile")
    print(img_hoechst.shape)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_mCherry, blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
    napari.run()

