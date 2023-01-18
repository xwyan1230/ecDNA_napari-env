import skimage.io as skio
import napari
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


sample = 'PC9'
# fov_lst = [fov_to_str(i) for i in np.arange(2, 37, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 7, 1)]  # Colo320DM
# fov_lst = [fov_to_str(i) for i in np.arange(1, 31, 1)]  # Colo320HSR
# fov_lst = ['MycI2_'+fov_to_str(i) for i in np.arange(1, 17, 1)]  # HCT116
fov_lst = list(np.arange(1, 9, 1)) + list(np.arange(10, 13, 1)) + [fov_to_str(i) for i in np.arange(13, 16, 1)] + [fov_to_str(i) for i in np.arange(17, 26, 1)] + [fov_to_str(i) for i in np.arange(27, 36, 1)] # PC3
# fov_lst = [fov_to_str(i) for i in np.arange(1, 2, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(1, 8, 1)] + ['Mn_'+fov_to_str(i) for i in np.arange(9, 19, 1)]  # PC9

for fov in range(len(fov_lst)):
    file_name = '%s_%s_Lng_SVCC_Processed001_RAW' % (sample, fov_lst[fov])
    img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
    img_RNAFISH = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_RNAFISH, blending='additive', colormap='red', contrast_limits=[0, img_RNAFISH.max()])
    # plt.imsave('%s%s/color_img/%s_%s_img_FISH_and_IF.tiff' % (output_dir, sample, sample, i), dis.blending(viewer))
    # viewer.close()
    napari.run()
