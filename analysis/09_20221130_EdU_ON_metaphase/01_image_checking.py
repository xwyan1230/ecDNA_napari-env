import skimage.io as skio
import napari
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221208_analysis_20221130_EdU_ON_metaphase/20221130_EdU_ON_metaphase/"
start_fov = 1
total_fov = 43


def fov_to_str(fov):
    if fov < 10:
        out = '00%s' % fov
    else:
        out = '0%s' % fov
    return out


# 2048x2048
for f in range(total_fov):
    fov = f + start_fov
    # file_name = 'Series%s_RAW' % fov_to_str(fov)
    file_name = '3144_RAW'
    img_hoechst = skio.imread("%s%s_ch00.tif" % (master_folder, file_name), plugin="tifffile")
    img_FISH = skio.imread("%s%s_ch02.tif" % (master_folder, file_name), plugin="tifffile")
    img_EdU = skio.imread("%s%s_ch01.tif" % (master_folder, file_name), plugin="tifffile")
    # img_EdU = np.concatenate([np.zeros(shape=[2048, 6]), img_EdU], axis=1) # for 2048x2048
    img_EdU = np.concatenate([np.zeros(shape=[3144, 8]), img_EdU], axis=1)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_FISH, blending='additive', colormap='green', contrast_limits=[0, img_FISH.max()])
    viewer.add_image(img_EdU, blending='additive', colormap='magenta', contrast_limits=[0, img_EdU.max()])
    napari.run()
