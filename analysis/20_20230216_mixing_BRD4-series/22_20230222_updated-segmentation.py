import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import matplotlib.pyplot as plt
import shared.display as dis
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.dataframe as dat
import tifffile as tif
from skimage import segmentation
import shared.objects as obj
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/singleZ/" % master_folder
data_dir2 = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

sample = 'DM-Ctrl_mix_mCh-Ctrl_forBRD3'
batch = 0
if batch > 0:
    file_name = '%s_%s_RAW' % (sample, batch)
else:
    file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")
# img_laminB_stack = skio.imread("%s%s/%s_ch03.tif" % (data_dir, sample, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 12
local_size = 200
rmax = 100
start_fov = 11


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out


for f in range(img_hoechst_stack.shape[0]):
    fov = start_fov + f
    print("Analyzing %s, fov %s" % (sample, fov))
    img_nuclear = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir2, sample, batch, sample, fov), plugin="tifffile")
    img_ecSeg = skio.imread("%s%s/seg_tif/%s/%s_%s_ecseg_n12.tif" % (data_dir2, sample, batch, sample, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))
    img_DNAFISH_seg1 = np.zeros_like(img_DNAFISH)

    # measure
    # get local images
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):

        print("Analyzing %s, fov %s, nuclear %s/%s" % (sample, fov + 1, i + 1, len(nuclear_props)))
        original_centroid_nuclear = nuclear_props[i].centroid
        position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, nuclear_props[i].label)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg = img_ecSeg.copy()
        local_DNAFISH_seg = local_DNAFISH_seg[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        local_nuclear_seg_mCherry = ima.img_local_seg(img_seg, position, nuclear_props[i].label)
        local_nuclear_props = regionprops(label(local_nuclear_seg))
        local_centroid = local_nuclear_props[0].centroid

        local_DNAFISH_singlet = local_DNAFISH.copy()
        local_DNAFISH_singlet[local_nuclear_seg == 0] = 0
        local_DNAFISH_seg1 = local_DNAFISH_singlet > 12000
        local_DNAFISH_seg1 = obj.remove_small(local_DNAFISH_seg1, 20)

        """viewer = napari.Viewer()
        # viewer.add_image(local_nuclear, blending='additive', colormap='blue', contrast_limits=[0, local_nuclear.max()])
        viewer.add_image(local_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, local_DNAFISH.max()])
        viewer.add_image(local_DNAFISH_seg, blending='additive', contrast_limits=[0, 1])
        viewer.add_image(local_DNAFISH_seg1, blending='additive', contrast_limits=[0, 1])
        napari.run()"""
        img_DNAFISH_seg1 = ima.image_paste_to(img_DNAFISH_seg1, local_DNAFISH_seg1,
                                             [int(original_centroid_nuclear[0] - local_centroid[0]),
                                              int(original_centroid_nuclear[1] - local_centroid[1])])

    tif.imwrite("%s%s/seg_tif/%s/%s_%s_ecseg1_n%s.tif" % (output_dir, sample, batch, sample, fov, n_nuclear_convex_dilation),
            img_DNAFISH_seg1)

    viewer = napari.Viewer()
    viewer.add_image(img_nuclear, blending='additive', colormap='blue', contrast_limits=[0, img_nuclear.max()])
    viewer.add_image(img_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, img_DNAFISH.max()])
    plt.imsave('%s%s/color_img/%s/%s_%s_img_womCh.tiff' % (output_dir, sample, batch, sample, fov), dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(dilation(img_seg, disk(2)), blending='additive', colormap='blue',
                     contrast_limits=[0, 1])
    viewer.add_image(img_DNAFISH_seg1, blending='additive', colormap='green', contrast_limits=[0, 1])
    plt.imsave('%s%s/color_img/%s/%s_%s_seg1_n%s.tiff' % (output_dir, sample, batch, sample, fov, n_nuclear_convex_dilation),
        dis.blending(viewer))
    viewer.close()

print("DONE!")