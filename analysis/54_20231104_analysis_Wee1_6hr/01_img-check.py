import skimage.io as skio
import napari
import os
import tifffile as tif
import shared.image as ima
from skimage.morphology import binary_dilation, disk
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231104_analysis_ColoDMandHSR_Wee1_6hr/"

sample = 'ColoHSR_Wee1_degrader_gH2AX_pH3'
conc = 'XY42'

# img
fov = 1  # 1, 3, 7, 9
img_CH1 = skio.imread("%s/data/%s/%s/Image_%s_CH1.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 2]
img_CH2 = skio.imread("%s/data/%s/%s/Image_%s_CH2.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 1]
img_CH3 = skio.imread("%s/data/%s/%s/Image_%s_CH3.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 0]
img_CH4 = skio.imread("%s/data/%s/%s/Image_%s_CH4.tif" % (master_folder, sample, conc, conc),
                      plugin="tifffile")[:, :, 0]
img_seg = np.zeros_like(img_CH1)

viewer = napari.Viewer()
viewer.add_image(img_CH1, blending='additive', colormap='blue', contrast_limits=[0, 30000])
viewer.add_image(img_CH2, blending='additive', colormap='green', contrast_limits=[0, 30000])
viewer.add_image(img_CH3, blending='additive', colormap='red', contrast_limits=[0, 30000])
viewer.add_image(img_CH4, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
shapes_add = viewer.add_shapes(name='add', ndim=2)
napari.run()

img_seg = ima.napari_add_or_remove_obj(shapes_add.data, 'add', img_seg)
img_seg = binary_dilation(img_seg, disk(7))

"""viewer = napari.Viewer()
viewer.add_image(img_CH1, blending='additive', colormap='blue', contrast_limits=[0, 30000])
viewer.add_image(img_CH2, blending='additive', colormap='green', contrast_limits=[0, 30000])
viewer.add_image(img_CH3, blending='additive', colormap='red', contrast_limits=[0, 30000])
viewer.add_image(img_CH4, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
viewer.add_image(img_seg, blending='additive', contrast_limits=[0, 1])
napari.run()"""

tif.imwrite("%s/data/%s/%s/Image_%s_seg.tif" % (master_folder, sample, conc, conc), img_seg)
