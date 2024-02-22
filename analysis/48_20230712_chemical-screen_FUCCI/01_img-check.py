import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_analysis_chemical-screen_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

plate = 'DM_FUCCI_24hr'
sample_lst = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C11', 'C10', 'C9', 'C8', 'C7', 'C6',
              'C5', 'C4', 'C3', 'C2', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'E11', 'E10',
              'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10',
              'F11', 'G11', 'G10', 'G9', 'G8', 'G7', 'G6', 'G5', 'G4', 'G3', 'G2']
sample = 'C3'
pos = sample_lst.index(sample)
print(pos)
for i in range(4):
    file_num = 4*pos+i+1
    if file_num < 10:
        file_folder = 'XY0%s' % file_num
    else:
        file_folder = 'XY%s' % file_num
    img_hoechst = skio.imread("%s%s/%s/Image_%s_CH1.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 2]  # 1536, 2048
    img_green = skio.imread("%s%s/%s/Image_%s_CH2.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 1]
    img_red = skio.imread("%s%s/%s/Image_%s_CH3.tif" % (data_dir, plate, file_folder, file_folder), plugin="tifffile")[..., 0]
    print(img_hoechst.shape)
    print(img_hoechst.max())
    # img_hoechst = np.concatenate([np.zeros(shape=[2048, 5]), img_hoechst], axis=1)[:2048, :2048]

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 255])
    viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, 255])
    viewer.add_image(img_red, blending='additive', colormap='red', contrast_limits=[0, 255])
    napari.run()
