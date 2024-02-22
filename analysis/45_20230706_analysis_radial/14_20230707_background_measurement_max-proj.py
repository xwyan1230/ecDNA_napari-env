import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.dataframe as dat
from skimage import segmentation
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/max_proj/" % master_folder
output_dir = "%sdata/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'
total_fov = 6
start_fov = 1

bg_lst = []

for f in range(total_fov):
    fov = start_fov + f
    img_DNAFISH_maxproj = skio.imread("%s%s/fov%s_DNAFISH_maxproj.tif" % (data_dir, sample, fov), plugin="tifffile")

    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH_maxproj, blending='additive', colormap='green', contrast_limits=[0, 65535])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()

    poly_data = shapes.data[0]
    shapes_layer = viewer.layers['Shapes']
    top, left = np.floor(np.min(poly_data, axis=0))
    bottom, right = np.ceil(np.max(poly_data, axis=0))
    top, bottom = np.clip((top, bottom), 0, img_DNAFISH_maxproj.shape[0] - 1).astype(int)
    left, right = np.clip((left, right), 0, img_DNAFISH_maxproj.shape[1] - 1).astype(int)
    output_shape = (bottom - top + 1, right - left + 1)
    # generate sub_masks and sub_channels
    sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[0]
    img_DNAFISH_mask = img_DNAFISH_maxproj[top:bottom + 1, left:right + 1] * sub_masks
    bg = regionprops(label(sub_masks), img_DNAFISH_mask)[0].intensity_mean
    bg_lst.append(bg)

print(bg_lst)
print(np.mean(bg_lst))

print("DONE!")