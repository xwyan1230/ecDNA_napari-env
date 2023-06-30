import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/singleZ/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-Ctrl_mix_mCh-BRD1'
batch = 0
if batch > 0:
    file_name = '%s_%s_RAW' % (sample, batch)
else:
    file_name = '%s_RAW' % sample
img_hoechst = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_DNAFISH = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_mCherry = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_laminB = skio.imread("%s%s/%s_ch03.tif" % (data_dir, sample, file_name), plugin="tifffile")

print(img_hoechst.shape)

if len(img_hoechst.shape) > 2:
    for i in range(img_hoechst.shape[0]):
        viewer = napari.Viewer()
        viewer.add_image(img_hoechst[i, :, :], blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
        viewer.add_image(img_DNAFISH[i, :, :], blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
        viewer.add_image(img_mCherry[i, :, :], blending='additive', colormap='red', contrast_limits=[0, img_mCherry.max()])
        viewer.add_image(img_laminB[i, :, :], blending='additive', colormap='magenta', contrast_limits=[0, img_laminB.max()])
        napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='red', contrast_limits=[0, 65535])
    viewer.add_image(img_mCherry, blending='additive', colormap='green', contrast_limits=[0, img_mCherry.max()])
    viewer.add_image(img_laminB, blending='additive', colormap='magenta', contrast_limits=[0, img_laminB.max()])
    napari.run()
