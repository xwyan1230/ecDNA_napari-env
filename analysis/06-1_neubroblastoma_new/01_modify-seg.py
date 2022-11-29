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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/E/"
sample = 'LY'
save_path = master_folder

# segmentation
local_factor_nuclear = 151
min_size_nuclear = 3000
max_size_nuclear = 15000
convex_conversion_threshold = 0.85

# other parameters
local_size = 100

# load images
img_nuclear = skio.imread("%s%s_b.BMP" % (master_folder, sample))[:, :, 2]
# 1024x1360
# img_nuclear = np.concatenate([img_nuclear, np.zeros(shape=[10, 1360])], axis=0)
# tif.imwrite("%s/%s_b_resize.tif" % (master_folder, sample), img_nuclear)
img_DNAFISH = skio.imread("%s%s_g.BMP" % (master_folder, sample))[:, :, 1]
# img_DNAFISH = np.concatenate([img_DNAFISH, np.zeros(shape=[10, 1360])], axis=0)
# tif.imwrite("%s/%s_g_resize.tif" % (master_folder, sample), img_DNAFISH)
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

# nuclear segmentation
img_nuclear_seg_convex = skio.imread("%s%s_seg.tif" % (master_folder, sample), plugin="tifffile")

# manual correction for nuclear segmentation
viewer = napari.Viewer()
viewer.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bg_corrected.max()])
viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green',
                 contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='viridis')
shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
napari.run()

img_nuclear_seg_convex = ima.napari_add_or_remove(shapes_seg_remove.data, 'remove', img_nuclear_seg_convex)
img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_seg_add.data, 'add', img_nuclear_seg_convex)
img_nuclear_seg_convex = obj.label_resort(img_nuclear_seg_convex)

img_nuclear_seg_convex[img_nuclear_seg_convex <= 11] = 0

viewer = napari.Viewer()
viewer.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear_bg_corrected.max()])
viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green',
                 contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='viridis')
shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
napari.run()

tif.imwrite("%s/%s_seg.tif" % (master_folder, sample), img_nuclear_seg_convex)
nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear_bg_corrected)

img_DNAFISH_seg = skio.imread("%s%s_ecSeg.tif" % (master_folder, sample), plugin="tifffile")

# manual correction for ecDNA segmentation
viewer = napari.Viewer()
viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green',
                 contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
viewer.add_image(img_DNAFISH_seg, blending='additive', contrast_limits=[0, 1])
shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
napari.run()

img_DNAFISH_seg = ima.napari_add_or_remove(shapes_seg_remove.data, 'remove', img_DNAFISH_seg)
img_DNAFISH_seg = ima.napari_add_or_remove(shapes_seg_add.data, 'add', img_DNAFISH_seg)

"""viewer = napari.Viewer()
viewer.add_image(local_nuclear, blending='additive', colormap='blue')
viewer.add_image(local_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
viewer.add_image(FISH_seg_local, blending='additive')
viewer.add_image(FISH_seg, blending='additive')
napari.run()"""

tif.imwrite("%s/%s_ecSeg.tif" % (master_folder, sample), img_DNAFISH_seg)

# viewer
viewer = napari.Viewer()
viewer.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
# viewer.add_image(img_centromere_bg_corrected, blending='additive', colormap='red')
viewer.add_image(img_nuclear_seg_convex, blending='additive')
plt.imsave('%s%s_nuclei.tiff' % (save_path, sample), dis.blending(viewer))
viewer.close()

viewer1 = napari.Viewer()
viewer1.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer1.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
viewer1.add_image(img_DNAFISH_seg, blending='additive')
plt.imsave('%s%s_DNAFISH.tiff' % (save_path, sample), dis.blending(viewer1))
viewer1.close()

viewer2 = napari.Viewer()
viewer2.add_image(img_nuclear_bg_corrected, blending='additive', colormap='blue')
viewer2.add_image(img_DNAFISH_bg_corrected, blending='additive', colormap='green')
plt.imsave('%s%s_img.tiff' % (save_path, sample), dis.blending(viewer2))
viewer2.close()

viewer3 = napari.Viewer()
viewer3.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue')
viewer3.add_image(img_DNAFISH_seg, blending='additive', colormap='green')
plt.imsave('%s%s_seg.tiff' % (save_path, sample), dis.blending(viewer3))
viewer3.close()

print("DONE!")
