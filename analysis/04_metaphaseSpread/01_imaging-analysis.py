import random
import shared.segmentation as seg
from skimage.filters import threshold_otsu
import shared.objects as obj
import tifffile as tif
import matplotlib.pyplot as plt
import skimage.io as skio
import napari
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220701_BRD4series/220630_Colo_Metaphasesa/"
sample = 'DM_Meta'

# IMAGING ANALYSIS
# load images
img_hoechst = skio.imread("%s%s/%s_Position 1_RAW_ch00.tif" % (master_folder, sample, sample), plugin="tifffile")
img_FISH = skio.imread("%s%s/%s_Position 1_RAW_ch01.tif" % (master_folder, sample, sample), plugin="tifffile")

viewer = napari.Viewer()
layer_image = viewer.add_image(img_FISH)
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]

top, left = np.floor(np.min(poly_data, axis=0))
bottom, right = np.ceil(np.max(poly_data, axis=0))
top, bottom = np.clip((top, bottom), 0, img_FISH.shape[0] - 1).astype(int)
left, right = np.clip((left, right), 0, img_FISH.shape[1] - 1).astype(int)
output_shape = (bottom - top + 1, right - left + 1)

shapes_layer = viewer.layers['Shapes']
sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[0]
cell_FISH = img_FISH[top:bottom+1, left:right+1] * sub_masks

viewer = napari.Viewer()
viewer.add_image(cell_FISH)
viewer.add_image(sub_masks)
napari.run()







