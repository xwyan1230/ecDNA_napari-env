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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220927_Jun_EGFR_RPAs33p_Edu/EGFR_RPAs33p_Edu/"
sample = 'gbm39ec hu'
master_path = '%s%s/' % (master_folder, sample)
end_fov = 10
start_fov = 6
total_fov = end_fov - start_fov + 1
save_path = master_path

for f in range(total_fov):
    fov = f + start_fov
    print("Analyzing %s, start nuclear segmentation FOV %s/%s" % (sample, fov, total_fov))
    file_prefix = "40x %s-%s" % (sample, fov)
    # LOAD IMAGE
    img_nuclear = skio.imread("%sC2-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_DNAFISH = skio.imread("%sC1-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_RPA = skio.imread("%sC3-%s.tif" % (master_path, file_prefix), plugin="tifffile")
    img_EdU = skio.imread("%sC4-%s.tif" % (master_path, file_prefix), plugin="tifffile")

    # nuclear segmentation
    img_nuclear_seg_convex = skio.imread('%s%s_nuclear_seg.tif' % (master_path, file_prefix), plugin="tifffile")

    # manual correction for nuclear segmentation
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='viridis')
    shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
    shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
    napari.run()

    img_nuclear_seg_convex = ima.napari_add_or_remove(shapes_seg_remove.data, 'remove', img_nuclear_seg_convex)
    img_nuclear_seg_convex = ima.napari_add_or_remove_obj(shapes_seg_add.data, 'add', img_nuclear_seg_convex)
    img_nuclear_seg_convex = obj.label_resort(img_nuclear_seg_convex)

    tif.imwrite('%s%s_nuclear_seg.tif' % (master_path, file_prefix), img_nuclear_seg_convex)

    # DNAFISH segmentation
    img_DNAFISH_seg = skio.imread('%s%s_DNAFISH_seg.tif' % (master_path, file_prefix), plugin="tifffile")

    # manual correction for ecDNA segmentation
    viewer = napari.Viewer()
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg, blending='additive', contrast_limits=[0, 1])
    shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
    shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
    napari.run()

    img_DNAFISH_seg = ima.napari_add_or_remove(shapes_seg_remove.data, 'remove', img_DNAFISH_seg)
    img_DNAFISH_seg = ima.napari_add_or_remove(shapes_seg_add.data, 'add', img_DNAFISH_seg)

    tif.imwrite('%s%s_DNAFISH_seg.tif' % (master_path, file_prefix), img_DNAFISH_seg)

    # RPA segmentation
    img_RPA_seg = skio.imread('%s%s_RPA_seg.tif' % (master_path, file_prefix), plugin="tifffile")

    # manual correction for ecDNA segmentation
    viewer = napari.Viewer()
    viewer.add_image(img_RPA, blending='additive', colormap='red', contrast_limits=[0, img_DNAFISH.max()])
    viewer.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue', contrast_limits=[0, 1])
    viewer.add_image(img_RPA_seg, blending='additive', contrast_limits=[0, 1])
    shapes_seg_remove = viewer.add_shapes(name='seg to be removed', ndim=2)
    shapes_seg_add = viewer.add_shapes(name='seg to be added', ndim=2)
    napari.run()

    img_RPA_seg = ima.napari_add_or_remove(shapes_seg_remove.data, 'remove', img_RPA_seg)
    img_RPA_seg = ima.napari_add_or_remove(shapes_seg_add.data, 'add', img_RPA_seg)

    tif.imwrite('%s%s_RPA_seg.tif' % (master_path, file_prefix), img_RPA_seg)

    # viewer
    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue')
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green')
    viewer.add_image(img_RPA, blending='additive', colormap='red')
    viewer.add_image(img_EdU, blending='additive', colormap='magenta')
    plt.imsave('%s%s_img.tiff' % (master_path, file_prefix), dis.blending(viewer))
    viewer.close()

    viewer1 = napari.Viewer()
    viewer1.add_image(img_nuclear_seg_convex, blending='additive', colormap='blue')
    viewer1.add_image(img_DNAFISH_seg, blending='additive', colormap='green', contrast_limits=[0, 1])
    viewer1.add_image(img_RPA_seg, blending='additive', colormap='red', contrast_limits=[0, 1])
    plt.imsave('%s%s_seg.tiff' % (master_path, file_prefix), dis.blending(viewer1))
    viewer1.close()

print("DONE!")
