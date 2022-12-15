import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-new'
file_name = 'DM_%s_acid_RAW' % sample
img_hoechst = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s%s_ch02.tif" % (data_dir, file_name), plugin="tifffile")
img_MYC = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")
# img_APC1 = skio.imread("%s%s_ch03.tif" % (master_folder, file_name), plugin="tifffile")
# img_EdU = np.concatenate([np.zeros(shape=[2048, 6]), img_EdU], axis=1)

print(img_hoechst.shape)

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, :, :], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH[i, :, :], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        viewer.add_image(img_MYC[i, :, :], blending='additive', colormap='red', contrast_limits=[0, img_MYC.max()])
        # viewer.add_image(img_APC1[i, :, :], blending='additive', colormap='magenta', contrast_limits=[0, img_APC1.max()])
        # plt.imsave('%s%s_%s_img_DNAFISH.tiff' % (output_dir, sample, i), dis.blending(viewer))
        napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_MYC, blending='additive', colormap='green', contrast_limits=[0, img_MYC.max()])
    # viewer.add_image(img_APC1, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    # plt.imsave('%s%s_img_647.tiff' % (output_dir, sample), dis.blending(viewer))
    napari.run()
