import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230508_analysis_Wee1-series_IF/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'Wee1_6hr'
img_hoechst = skio.imread("%s/%s/Hoechst.tif" % (data_dir, sample), plugin="tifffile")
img_NPC = skio.imread("%s/%s/NPC.tif" % (data_dir, sample), plugin="tifffile")
img_laminB1 = skio.imread("%s/%s/laminB1.tif" % (data_dir, sample), plugin="tifffile")

print(img_hoechst.shape)

for i in range(img_hoechst.shape[0]):
    viewer = napari.Viewer()
    viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
    viewer.add_image(img_NPC, blending='additive', colormap='red', contrast_limits=[0, img_NPC.max()])
    viewer.add_image(img_laminB1, blending='additive', colormap='magenta', contrast_limits=[0, img_laminB1.max()])
    shapes = viewer.add_shapes()
    napari.run()

