import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import shared.objects as obj
import shared.image as ima
from skimage import segmentation
import numpy as np
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import pandas as pd
import tifffile as tif
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/A/"
sample = 'EB'
save_path = master_folder

"""# load images
img_nuclear_seg_binary = skio.imread("%s%s_seg.tif" % (master_folder, sample), plugin="tifffile")
img_nuclear_seg_convex = label(img_nuclear_seg_binary)
tif.imwrite("%s/%s_seg.tif" % (master_folder, sample), img_nuclear_seg_convex)"""

# load images
img_nuclear = skio.imread("%s%s_b.BMP" % (master_folder, sample))[:, :, 2]
img_DNAFISH = skio.imread("%s%s_g.BMP" % (master_folder, sample))[:, :, 1]
img_centromere = skio.imread("%s%s_r.BMP" % (master_folder, sample))[:, :, 0]

# bg_correction
bg_val_nuclear = seg.get_bg_int([img_nuclear])[0]
bg_val_DNAFISH = seg.get_bg_int([img_DNAFISH])[0]
bg_val_centromere = seg.get_bg_int([img_centromere])[0]
img_nuclear_bg_corrected = img_nuclear.astype(float) - np.ones_like(img_nuclear) * bg_val_nuclear
img_nuclear_bg_corrected[img_nuclear_bg_corrected < 0] = 0
img_DNAFISH_bg_corrected = img_DNAFISH.astype(float) - np.ones_like(img_DNAFISH) * bg_val_DNAFISH
img_DNAFISH_bg_corrected[img_DNAFISH_bg_corrected < 0] = 0
img_centromere_bg_corrected = img_centromere.astype(float) - np.ones_like(img_centromere) * bg_val_centromere
img_centromere_bg_corrected[img_centromere_bg_corrected < 0] = 0

img_seg = skio.imread("%s%s_seg.tif" % (master_folder, sample), plugin="tifffile")
img_seg_updated = img_seg
ecSeg = skio.imread("%s%s_ecSeg.tif" % (master_folder, sample), plugin="tifffile")
ecSeg_updated = binary_erosion(ecSeg)
tif.imwrite("%s/%s_ecSeg.tif" % (master_folder, sample), ecSeg_updated)

# viewer
viewer = napari.Viewer()
viewer.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
# viewer.add_image(img_centromere_bg_corrected, blending='additive', colormap='red')
viewer.add_image(img_seg_updated, blending='additive')
plt.imsave('%s%s_nuclei.tiff' % (save_path, sample), dis.blending(viewer))
viewer.close()

viewer1 = napari.Viewer()
viewer1.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer1.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
viewer1.add_image(ecSeg_updated, blending='additive')
plt.imsave('%s%s_DNAFISH.tiff' % (save_path, sample), dis.blending(viewer1))
viewer1.close()

viewer2 = napari.Viewer()
viewer2.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer2.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
plt.imsave('%s%s_img.tiff' % (save_path, sample), dis.blending(viewer2))
viewer2.close()

viewer3 = napari.Viewer()
viewer3.add_image(img_seg_updated, blending='additive', colormap='blue')
viewer3.add_image(ecSeg_updated, blending='additive', colormap='green')
plt.imsave('%s%s_seg.tiff' % (save_path, sample), dis.blending(viewer3))
viewer3.close()

print("DONE!")