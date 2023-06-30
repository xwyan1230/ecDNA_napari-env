import skimage.io as skio
import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import math
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
from skimage.measure import label, regionprops_table, regionprops
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'





img_hoechst = skio.imread("%s%s/DM_hoechst_IF_22000-33601.tif" % (data_dir, sample), plugin="tifffile")
img_nuclear_seg = skio.imread("%s%s/DM_nuclear_red_green_IF_22000-33601_seg_convex.tif" % (data_dir, sample), plugin="tifffile")

props = regionprops(img_nuclear_seg, img_hoechst)

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_nuclear_seg, blending='additive', contrast_limits=[0, 1])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
shapes_layer = viewer.layers['Shapes']

label_lst = []
for i in range(len(shapes.data)):
    poly_data = shapes.data[i]

    # ***** generate local images *****
    # reshape sub_masks
    top, left = np.floor(np.min(poly_data, axis=0))
    bottom, right = np.ceil(np.max(poly_data, axis=0))
    top, bottom = np.clip((top, bottom), 0, img_hoechst.shape[0] - 1).astype(int)
    left, right = np.clip((left, right), 0, img_hoechst.shape[1] - 1).astype(int)
    output_shape = (bottom - top + 1, right - left + 1)
    # generate sub_masks and sub_channels
    sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[i]
    nuclear_seg_mask = img_nuclear_seg[top:bottom + 1, left:right + 1] * sub_masks
    print(nuclear_seg_mask.max())
    label_lst.append(nuclear_seg_mask.max())

data = pd.DataFrame(columns=['label'])
data['label'] = label_lst

data.to_csv('%s%s_mitosis_label.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")



