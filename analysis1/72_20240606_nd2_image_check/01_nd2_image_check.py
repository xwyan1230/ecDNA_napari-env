import nd2
import napari
import pandas as pd
import numpy as np
import random
import shared.image as ima
import tifffile as tif
import skimage.io as skio
import matplotlib.pyplot as plt
import shared.display as dis

print(np.arange(1, 61, 1))
print(random.sample(list(np.arange(1, 61, 1)), 8))

"""# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240606_nd2/"

sample = 'VT40_002'
img_stack = nd2.imread('%s%s.nd2' % (master_folder, sample))
print(img_stack.shape)

# img = skio.imread("%sImage_XY01_00001_CH2.tif" % master_folder, plugin="tifffile")
# print(img.shape)

img1_stack = img_stack[:, :, :]
# img1_stack = img_stack[:, 0, :, :]
# img2_stack = img_stack[:, 1, :, :]

img1 = img_stack[5, :, :]
# img1 = img_stack[3, 0, :, :]
# img2 = img_stack[3, 1, :, :]
# img3 = img_stack[5, 2, :, :]

viewer = napari.Viewer()
viewer.add_image(img1_stack, blending='additive', contrast_limits=[0, 25000])
# viewer.add_image(img1_stack, blending='additive', colormap='blue', contrast_limits=[0, 5000])
# viewer.add_image(img2_stack, blending='additive', colormap='green', contrast_limits=[0, 60000])
# viewer.add_image(img3, blending='additive', colormap='red', contrast_limits=[0, 5000])
# plt.imsave("%s%s.tiff" % (master_folder, sample), dis.blending(viewer))
napari.run()

viewer = napari.Viewer()
viewer.add_image(img1, blending='additive', contrast_limits=[0, 23000])
# viewer.add_image(img1, blending='additive', colormap='blue', contrast_limits=[0, 5000])
# viewer.add_image(img2, blending='additive', colormap='green', contrast_limits=[0, 60000])
# viewer.add_image(img3, blending='additive', colormap='red', contrast_limits=[0, 5000])
plt.imsave("%s%s.tiff" % (master_folder, sample), dis.blending(viewer))
viewer.close()
"""

"""import matplotlib
k = 200
cmap = matplotlib.cm.get_cmap('viridis')

rgba = cmap(0.5)[:3]
print(rgba)"""