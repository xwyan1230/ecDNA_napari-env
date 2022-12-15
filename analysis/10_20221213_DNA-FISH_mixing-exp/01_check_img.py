import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H2B+BRD4'
file_name = '%s_RAW' % sample
img_hoechst = skio.imread("%s%s_ch00.tif" % (data_dir, file_name), plugin="tifffile")
img_PE = skio.imread("%s%s_ch01.tif" % (data_dir, file_name), plugin="tifffile")
img_APC = skio.imread("%s%s_ch02.tif" % (data_dir, file_name), plugin="tifffile")
# img_APC1 = skio.imread("%s%s_ch03.tif" % (master_folder, file_name), plugin="tifffile")
# img_EdU = np.concatenate([np.zeros(shape=[2048, 6]), img_EdU], axis=1)

print(img_hoechst.shape)

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, :, :], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_PE[i, :, :], blending='additive', colormap='green', contrast_limits=[0, img_PE.max()])
        viewer.add_image(img_APC[i, :, :], blending='additive', colormap='red', contrast_limits=[0, img_APC.max()])
        # viewer.add_image(img_APC1[i, :, :], blending='additive', colormap='magenta', contrast_limits=[0, img_APC1.max()])
        # plt.imsave('%s%s_%s_img_DNAFISH.tiff' % (output_dir, sample, i), dis.blending(viewer))
        napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_PE, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_APC, blending='additive', colormap='green', contrast_limits=[0, img_APC.max()])
    # viewer.add_image(img_APC1, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
    # plt.imsave('%s%s_img_647.tiff' % (output_dir, sample), dis.blending(viewer))
    napari.run()
