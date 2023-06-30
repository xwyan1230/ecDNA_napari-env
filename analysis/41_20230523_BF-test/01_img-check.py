import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, threshold_yen, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import shared.objects as obj
import shared.segmentation as seg
from scipy import ndimage
import shared.image as ima

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230524_analysis_BF-test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'test'
img_BF = ima.img_to_int(skio.imread("%s/%s_CH4.tif" % (data_dir, sample), plugin="tifffile"))
img_green = skio.imread("%s/%s_CH2.tif" % (data_dir, sample), plugin="tifffile")[:, :, 1]
img_red = skio.imread("%s/%s_CH3.tif" % (data_dir, sample), plugin="tifffile")[:, :, 0]
img_hoechst = skio.imread("%s/%s_CH1.tif" % (data_dir, sample), plugin="tifffile")[:, :, 2]

local = threshold_local(img_green, 99)
img_seg = img_green > local
img_seg = binary_erosion(img_seg, disk(2))

"""local = threshold_local(img_BF, 99)
img_BF_seg = img_BF > local
img_BF_seg = binary_erosion(img_BF_seg)

img_BF_seg1 = binary_dilation(img_BF_seg)
img_BF_seg_final1 = ndimage.binary_fill_holes(img_BF_seg1)
img_BF_seg_final1[img_BF_seg1 == 1] = 0
img_BF_seg_final1 = obj.remove_small(img_BF_seg_final1, min_size=150)
img_BF_seg_final1 = obj.remove_large(img_BF_seg_final1, max_size=700)

img_BF_seg2 = binary_dilation(img_BF_seg, disk(2))
img_BF_seg_final2 = ndimage.binary_fill_holes(img_BF_seg2)
img_BF_seg_final2[img_BF_seg2 == 1] = 0
img_BF_seg_final2 = obj.remove_small(img_BF_seg_final2, min_size=150)
img_BF_seg_final2 = obj.remove_large(img_BF_seg_final2, max_size=700)"""

viewer = napari.Viewer()
viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, img_hoechst.max()])
viewer.add_image(img_green, blending='additive', colormap='green', contrast_limits=[0, img_green.max()])
viewer.add_image(img_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
"""viewer.add_image(img_BF, blending='additive', contrast_limits=[0, img_BF.max()])
# viewer.add_image(img_BF_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
viewer.add_image(img_BF_seg_final1, blending='additive', colormap='blue', contrast_limits=[0, 1])
viewer.add_image(img_BF_seg_final2, blending='additive', colormap='blue', contrast_limits=[0, 1])"""
napari.run()

