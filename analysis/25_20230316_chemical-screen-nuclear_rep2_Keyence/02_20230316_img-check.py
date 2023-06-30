import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230316_analysis_chemical-screen-nuclear_rep2_Keyence/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

plate = 'HSR_24hr'
num_total_samples = 240
start_num = 0
for f in range(num_total_samples):
    fov = start_num + f
    file_name = 'Image_XY0%s' % (fov+1) if fov<9 else 'Image_XY%s' % (fov+1)
    print(file_name)
    img_hoechst = skio.imread("%s%s/raw_img/%s_CH1.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 2] # 1440x1920
    img_MYC = skio.imread("%s%s/raw_img/%s_CH3.tif" % (data_dir, plate, file_name), plugin="tifffile")[:, :, 0]
    # img_hoechst = np.concatenate([np.zeros(shape=[2048, 5]), img_hoechst], axis=1)[:2048, :2048]
    print(img_MYC.shape)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
    viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 255])
    napari.run()
