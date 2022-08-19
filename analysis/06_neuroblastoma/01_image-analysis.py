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
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220817_neuroblastoma/"
sample = 'NG'
save_path = master_folder

# segmentation
local_factor_nuclear = 151
min_size_nuclear = 3000
max_size_nuclear = 15000
convex_conversion_threshold = 0.85

# other parameters
local_size = 100

# load images
img_nuclear = skio.imread("%s%s_2b.BMP" % (master_folder, sample))[:, :, 2]
img_DNAFISH = skio.imread("%s%s_2g.BMP" % (master_folder, sample))[:, :, 1]
img_centromere = skio.imread("%s%s_2r.BMP" % (master_folder, sample))[:, :, 0]

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
img_nuclear_seg = seg.nuclear_seg(img_nuclear_bg_corrected, local_factor=local_factor_nuclear,
                                  min_size=min_size_nuclear, max_size=max_size_nuclear)
img_nuclear_seg_convex = obj.label_resort(seg.obj_to_convex_filter(img_nuclear_seg,
                                                                   threshold=convex_conversion_threshold))
nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear_bg_corrected)

img_DNAFISH_seg = np.zeros_like(img_DNAFISH_bg_corrected)

for i in range(len(nuclear_props)):
    original_centroid_nuclear = nuclear_props[i].centroid
    position = ima.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = ima.img_local_seg(img_nuclear_seg_convex, position, nuclear_props[i].label)
    local_nuclear = img_nuclear_bg_corrected.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH_bg_corrected.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]

    # ecDNA segmentation
    local_DNAFISH_singlet = local_DNAFISH.copy()
    local_DNAFISH_singlet[local_nuclear_seg_convex == 0] = 0
    otsu_threshold_val_local_DNAFISH = threshold_otsu(local_DNAFISH_singlet)

    threshold_min = otsu_threshold_val_local_DNAFISH + (img_DNAFISH_bg_corrected.max() - otsu_threshold_val_local_DNAFISH) / 2
    threshold_min7 = otsu_threshold_val_local_DNAFISH + (img_DNAFISH_bg_corrected.max() - otsu_threshold_val_local_DNAFISH) / 5

    FISH_seg_local = np.zeros_like(local_DNAFISH)

    for k in range(10):
        local = threshold_local(local_DNAFISH, 10*k+7)
        out = local_DNAFISH_singlet > local
        out = binary_erosion(out)
        if k == 0:
            out = binary_dilation(out)
            out_label = label(out)
            out_props = regionprops(out_label, local_DNAFISH_singlet)
            for j in range(len(out_props)):
                temp = np.zeros_like(local_DNAFISH)
                temp[out_label == out_props[j].label] = 1
                temp_outer_edge = binary_dilation(temp, disk(4))
                temp_outer_edge[temp == 1] = 0
                mean_int_outer_edge = np.sum(local_DNAFISH * temp_outer_edge)/np.sum(temp_outer_edge)
                if (out_props[j].intensity_mean/mean_int_outer_edge > 1.2) & (out_props[j].area > 5) & \
                        (out_props[j].intensity_mean > threshold_min7):
                    FISH_seg_local[out_label == out_props[j].label] = 1
        else:
            out_label = label(out)
            out_props = regionprops(out_label, local_DNAFISH_singlet)
            for j in range(len(out_props)):
                if out_props[j].intensity_mean > threshold_min:
                    FISH_seg_local[out_label == out_props[j].label] = 1

    bg_val = otsu_threshold_val_local_DNAFISH * 3
    extreme_val = local_DNAFISH_singlet.max() * 2 / otsu_threshold_val_local_DNAFISH
    maxima = extrema.h_maxima(local_DNAFISH, extreme_val)
    elevation_map = sobel(local_DNAFISH)
    markers = np.zeros_like(local_DNAFISH)
    markers[local_DNAFISH_singlet < bg_val] = 1
    markers[maxima == 1] = 2
    seg_wat = segmentation.watershed(elevation_map, markers)
    seg_wat_filter = obj.label_remove_small(label(seg_wat), min_size=6)
    FISH_seg_watershed = np.zeros_like(seg_wat)
    FISH_seg_watershed[seg_wat_filter > 1] = 1

    FISH_seg = FISH_seg_watershed.copy()
    FISH_seg[FISH_seg_local == 1] = 1
    FISH_seg[local_nuclear_seg_convex == 0] = 0

    img_DNAFISH_seg = ima.image_paste_to(img_DNAFISH_seg, FISH_seg, [int(original_centroid_nuclear[0]-100), int(original_centroid_nuclear[1]-100)])

    """viewer = napari.Viewer()
    viewer.add_image(local_nuclear, blending='additive', colormap='blue')
    viewer.add_image(local_DNAFISH, blending='additive', colormap='green',
                     contrast_limits=[0, img_DNAFISH_bg_corrected.max()])
    viewer.add_image(FISH_seg_local, blending='additive')
    viewer.add_image(FISH_seg, blending='additive')
    napari.run()"""

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
