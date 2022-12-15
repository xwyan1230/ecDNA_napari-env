import skimage.io as skio
import pandas as pd
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import tifffile as tif
import shared.dataframe as dat
import matplotlib.pyplot as plt
import shared.display as dis
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir1 = "%sdata/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

mCherry_pos_color = 'gray'
mCherry_neg_color = 'red'
mCherry_pos_threshold = 10000
mCherry_neg_threshold = 5000

sample = 'DM_mix_DM-H2B-mCherry'
file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s_ch00.tif" % (data_dir1, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s_ch01.tif" % (data_dir1, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s_ch02.tif" % (data_dir1, file_name), plugin="tifffile")

df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir2, sample)), na_values=['.'], sep='\t')

n_nuclear_convex_dilation = 2
local_size = 200
size = 20*local_size + 22


def grid_location(num: int):
    row = int(num / 10)
    column = num % 10
    x = row * (2 * local_size + 2) + 2 + local_size
    y = column * (2 * local_size + 2) + 2 + local_size
    return x, y


img_pos_count = 0
img_neg_count = 0
pos_count = 0
neg_count = 0
img_mCherry_pos_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_pos_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

for fov in range(img_mCherry_stack.shape[0]):
    if pos_count == 100:
        if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
            os.makedirs("%s%s/grid/" % (output_dir, sample))
        tif.imwrite("%s%s/grid/%s_pos_hoechst_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_hoechst)
        tif.imwrite("%s%s/grid/%s_pos_DNAFISH_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_DNAFISH)

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        plt.imsave('%s%s/grid/%s_pos_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        plt.imsave('%s%s/grid/%s_pos_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        plt.imsave('%s%s/grid/%s_pos_hoechst_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        pos_count = 0
        img_pos_count += 1
        img_mCherry_pos_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_pos_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

    if neg_count == 100:
        if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
            os.makedirs("%s%s/grid/" % (output_dir, sample))
        tif.imwrite("%s%s/grid/%s_neg_hoechst_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_hoechst)
        tif.imwrite("%s%s/grid/%s_neg_DNAFISH_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_DNAFISH)

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        plt.imsave('%s%s/grid/%s_neg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green',
                         contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        plt.imsave('%s%s/grid/%s_neg_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue',
                         contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        plt.imsave('%s%s/grid/%s_neg_hoechst_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        neg_count = 0
        img_neg_count += 1
        img_mCherry_neg_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_neg_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)

    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s_%s_seg.tif" % (data_dir2, sample, sample, fov), plugin="tifffile")

    mCherry_props = regionprops(img_seg, img_mCherry)

    for i in range(len(mCherry_props)):
        original_centroid_nuclear = mCherry_props[i].centroid
        position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_seg, position, mCherry_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_centroid = local_nuclear_props[0].centroid
        local_seg_filter = dilation(local_nuclear_seg_convex, disk(10))
        local_nuclear[local_seg_filter == 0] = 0
        local_DNAFISH[local_seg_filter == 0] = 0

        if mCherry_props[i].intensity_mean > mCherry_pos_threshold:
            if pos_count < 100:
                paste_to_centroid = [grid_location(pos_count)[0], grid_location(pos_count)[1]]
                print(paste_to_centroid)
                pos_count += 1
                img_mCherry_pos_hoechst = ima.image_paste_to(img_mCherry_pos_hoechst, local_nuclear,
                                                 [int(paste_to_centroid[0] - local_centroid[0]),
                                                  int(paste_to_centroid[1] - local_centroid[1])])
                img_mCherry_pos_DNAFISH = ima.image_paste_to(img_mCherry_pos_DNAFISH, local_DNAFISH,
                                                             [int(paste_to_centroid[0] - local_centroid[0]),
                                                              int(paste_to_centroid[1] - local_centroid[1])])
        elif mCherry_props[i].intensity_mean < mCherry_neg_threshold:
            if neg_count < 100:
                paste_to_centroid = [grid_location(neg_count)[0], grid_location(neg_count)[1]]
                print(paste_to_centroid)
                neg_count += 1
                img_mCherry_neg_hoechst = ima.image_paste_to(img_mCherry_neg_hoechst, local_nuclear,
                                                             [int(paste_to_centroid[0] - local_centroid[0]),
                                                              int(paste_to_centroid[1] - local_centroid[1])])
                img_mCherry_neg_DNAFISH = ima.image_paste_to(img_mCherry_neg_DNAFISH, local_DNAFISH,
                                                             [int(paste_to_centroid[0] - local_centroid[0]),
                                                              int(paste_to_centroid[1] - local_centroid[1])])

        """viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        # plt.imsave('%s%s/color_img/%s_%s_img_label.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
        napari.run()

        viewer1 = napari.Viewer()
        viewer1.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
        viewer1.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
        # plt.imsave('%s%s/color_img/%s_%s_img_label.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
        napari.run()"""

if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
    os.makedirs("%s%s/grid/" % (output_dir, sample))
tif.imwrite("%s%s/grid/%s_pos_hoechst_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_hoechst)
tif.imwrite("%s%s/grid/%s_pos_DNAFISH_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_DNAFISH)
tif.imwrite("%s%s/grid/%s_neg_hoechst_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_hoechst)
tif.imwrite("%s%s/grid/%s_neg_DNAFISH_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_DNAFISH)

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
plt.imsave('%s%s/grid/%s_pos_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
plt.imsave('%s%s/grid/%s_pos_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
plt.imsave('%s%s/grid/%s_pos_hoechst_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
plt.imsave('%s%s/grid/%s_neg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green',
                 contrast_limits=[0, max(img_mCherry_pos_DNAFISH.max(), img_mCherry_neg_DNAFISH.max())])
plt.imsave('%s%s/grid/%s_neg_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue',
                 contrast_limits=[0, max(img_mCherry_pos_hoechst.max(), img_mCherry_neg_hoechst.max())])
plt.imsave('%s%s/grid/%s_neg_hoechst_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()