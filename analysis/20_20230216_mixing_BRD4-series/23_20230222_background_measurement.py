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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/singleZ/" % master_folder
data_dir2 = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

sample1 = 'DM-Ctrl_mix_mCh-Ctrl'
figure_name = sample1
pos_threshold = 15000
neg_threshold = 12000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'
new_seg = 11000

start_fov = 15

batch = 2
if batch > 0:
    file_name = '%s_%s_RAW' % (sample1, batch)
else:
    file_name = '%s_RAW' % sample1
img_hoechst_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample1, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample1, file_name), plugin="tifffile")
# img_laminB_stack = skio.imread("%s%s/%s_ch03.tif" % (data_dir, sample, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 12
local_size = 200
rmax = 100
max_area = 60000
min_circ = 0.8


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask.flatten()) if e != 0]
    out = list(map(img.flatten().__getitem__, index))
    return out

pos_bg = []
neg_bg = []

for f in range(img_hoechst_stack.shape[0]):
    fov = start_fov+f
    print("Analyzing %s, fov %s" % (sample1, fov))
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry_bgc = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir2, sample1, batch, sample1, fov), plugin="tifffile")

    if n_nuclear_convex_dilation > 0:
        img_nuclear_seg = dilation(img_seg, disk(n_nuclear_convex_dilation))

    # measure
    # get local images
    mCherry_props = regionprops(img_seg, img_mCherry_bgc)
    for i in range(min(len(mCherry_props), 5)):
        if (mCherry_props[i].area < max_area) & (
                (4 * math.pi * mCherry_props[i].area) / (mCherry_props[i].perimeter ** 2) > min_circ):

            print("Analyzing %s, fov %s, nuclear %s/%s" % (sample1, fov, i + 1, len(mCherry_props)))
            original_centroid_nuclear = mCherry_props[i].centroid
            position = ima.img_local_position(img_nuclear_seg, original_centroid_nuclear, local_size)
            local_nuclear_seg = ima.img_local_seg(img_nuclear_seg, position, mCherry_props[i].label)
            local_nuclear = img_nuclear_bgc.copy()
            local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH = img_DNAFISH_bgc.copy()
            local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
            local_DNAFISH[local_nuclear_seg == 0] = 0

            viewer = napari.Viewer()
            viewer.add_image(local_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, local_DNAFISH.max()])
            shapes = viewer.add_shapes(name='Shapes', ndim=2)
            napari.run()

            poly_data = shapes.data[0]
            shapes_layer = viewer.layers['Shapes']
            top, left = np.floor(np.min(poly_data, axis=0))
            bottom, right = np.ceil(np.max(poly_data, axis=0))
            top, bottom = np.clip((top, bottom), 0, local_DNAFISH.shape[0] - 1).astype(int)
            left, right = np.clip((left, right), 0, local_DNAFISH.shape[1] - 1).astype(int)
            output_shape = (bottom - top + 1, right - left + 1)
            # generate sub_masks and sub_channels
            sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[0]
            local_DNAFISH_mask = local_DNAFISH[top:bottom + 1, left:right + 1] * sub_masks
            bg = regionprops(label(sub_masks), local_DNAFISH_mask)[0].intensity_mean

            if mCherry_props[i].intensity_mean > pos_threshold:
                pos_bg.append(bg)
                print('pos: %s' % len(pos_bg))
            elif mCherry_props[i].intensity_mean < neg_threshold:
                neg_bg.append(bg)
                print('neg: %s' % len(neg_bg))

    background = pd.DataFrame()
    background['sample'] = [sample1_pos] * len(pos_bg) + [sample1_neg] * len(neg_bg)
    background['bg'] = pos_bg + neg_bg
    if batch > 0:
        background.to_csv('%s%s_%s_background.txt' % (output_dir, sample1, batch), index=False, sep='\t')
    else:
        background.to_csv('%s%s_background.txt' % (output_dir, sample1), index=False, sep='\t')


print("DONE!")