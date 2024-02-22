import skimage.io as skio
import napari
import imutils
import shared.image as ima
import numpy as np
import tifffile as tif
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_Natasha_FFPE_tissue/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

prefix = '108_D2_GRmyc_REDccnd1 Merged_Lng_LVCC_Processed001_RAW_'
# prefix = '110_C3_GRmyc_REDccnd1 Merged_Lng_LVCC_Processed001_RAW_'

img_3d_nuclear = skio.imread("%s%sch00.tif" % (data_dir, prefix), plugin="tifffile")
img_3d_DNAFISH = skio.imread("%s%sch01.tif" % (data_dir, prefix), plugin="tifffile")
img_3d_laminB = skio.imread("%s%sch02.tif" % (data_dir, prefix), plugin="tifffile")

viewer = napari.Viewer()
viewer.add_image(img_3d_nuclear, blending='additive', colormap='blue', name='nuclei')
viewer.add_image(img_3d_DNAFISH, blending='additive', colormap='green', name='DNAFISH')
viewer.add_image(img_3d_laminB, blending='additive', colormap='red', name='laminB')
# viewer.window.add_plugin_dock_widget(plugin_name='napari-animation')
napari.run()

