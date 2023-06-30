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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sdata/singleZ/" % master_folder
data_dir2 = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH1/" % master_folder

mCherry_pos_color = 'gray'
mCherry_neg_color = 'red'
mCherry_pos_threshold = 21000
mCherry_neg_threshold = 16000

sample = 'DM-Ctrl_mix_mCh-BRD4'
batch = 0
if batch > 0:
    file_name = '%s_%s_RAW' % (sample, batch)
else:
    file_name = '%s_RAW' % sample
img_hoechst_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%s%s/%s_ch02.tif" % (data_dir, sample, file_name), plugin="tifffile")
img_mCherry_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, sample, file_name), plugin="tifffile")

n_nuclear_convex_dilation = 4
df = pd.read_csv(("%s%s_n%s.txt" % (data_dir2, sample, n_nuclear_convex_dilation)), na_values=['.'], sep='\t')

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
img_mCherry_pos_mCherry = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_pos_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_pos_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_mCherry = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_pos_seg = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_pos_ecseg = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_seg = np.zeros(shape=(size, size), dtype=np.uint16)
img_mCherry_neg_ecseg = np.zeros(shape=(size, size), dtype=np.uint16)

for fov in range(img_mCherry_stack.shape[0]):
    if pos_count == 100:
        if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
            os.makedirs("%s%s/grid/" % (output_dir, sample))
        tif.imwrite("%s%s/grid/%s_pos_hoechst_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_hoechst)
        tif.imwrite("%s%s/grid/%s_pos_DNAFISH_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_DNAFISH)
        tif.imwrite("%s%s/grid/%s_pos_seg_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_seg)
        tif.imwrite("%s%s/grid/%s_pos_ecseg_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_ecseg)
        tif.imwrite("%s%s/grid/%s_pos_mCherry_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_mCherry)

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_pos_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_pos_withmCh_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_pos_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
        viewer.add_image(img_mCherry_pos_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
        plt.imsave('%s%s/grid/%s_pos_seg_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_pos_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_pos_hoechst_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
        viewer.close()

        pos_count = 0
        img_pos_count += 1
        img_mCherry_pos_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_pos_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_pos_seg = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_pos_ecseg = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_pos_mCherry = np.zeros(shape=(size, size), dtype=np.uint16)

    if neg_count == 100:
        if not os.path.exists("%s%s/grid/" % (output_dir, sample)):
            os.makedirs("%s%s/grid/" % (output_dir, sample))
        tif.imwrite("%s%s/grid/%s_neg_hoechst_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_hoechst)
        tif.imwrite("%s%s/grid/%s_neg_DNAFISH_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_DNAFISH)
        tif.imwrite("%s%s/grid/%s_neg_seg_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_seg)
        tif.imwrite("%s%s/grid/%s_neg_ecseg_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_ecseg)
        tif.imwrite("%s%s/grid/%s_neg_mCherry_%s.tif" % (output_dir, sample, sample, img_neg_count),
                    img_mCherry_neg_mCherry)

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_neg_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_neg_withmCh_%s.tiff' % (output_dir, sample, sample, img_neg_count),
                   dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_neg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
        viewer.add_image(img_mCherry_neg_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
        plt.imsave('%s%s/grid/%s_neg_seg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_neg_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        viewer = napari.Viewer()
        viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        plt.imsave('%s%s/grid/%s_neg_hoechst_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
        viewer.close()

        neg_count = 0
        img_neg_count += 1
        img_mCherry_neg_hoechst = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_neg_DNAFISH = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_neg_seg = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_neg_ecseg = np.zeros(shape=(size, size), dtype=np.uint16)
        img_mCherry_neg_mCherry = np.zeros(shape=(size, size), dtype=np.uint16)

    print(fov)
    img_nuclear_bgc = img_hoechst_stack[fov, :, :]
    img_mCherry = img_mCherry_stack[fov, :, :]
    img_DNAFISH_bgc = img_DNAFISH_stack[fov, :, :]
    img_seg = skio.imread("%s%s/seg_tif/%s/%s_%s_seg.tif" % (data_dir2, sample, batch, sample, fov), plugin="tifffile")
    img_ecseg = skio.imread("%s%s/seg_tif/%s/%s_%s_ecseg_n12.tif" % (data_dir2, sample, batch, sample, fov),
                            plugin="tifffile")

    mCherry_props = regionprops(img_seg, img_mCherry)

    for i in range(len(mCherry_props)):
        original_centroid_nuclear = mCherry_props[i].centroid
        position = ima.img_local_position(img_seg, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = ima.img_local_seg(img_seg, position, mCherry_props[i].label)
        local_nuclear = img_nuclear_bgc.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH_bgc.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_mCherry = img_mCherry.copy()
        local_mCherry = local_mCherry[position[0]:position[1], position[2]:position[3]]
        local_ecseg = img_ecseg.copy()
        local_ecseg = local_ecseg[position[0]:position[1], position[2]:position[3]]
        local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
        local_DNAFISH_props = regionprops(label(local_nuclear_seg_convex), local_DNAFISH)
        local_centroid = local_nuclear_props[0].centroid
        local_seg_filter = dilation(local_nuclear_seg_convex, disk(10))
        local_nuclear[local_seg_filter == 0] = 0
        local_DNAFISH[local_seg_filter == 0] = 0
        local_ecseg[local_seg_filter == 0] = 0
        local_mCherry[local_seg_filter == 0] = 0

        total_int_DNAFISH = local_nuclear_props[0].area * local_DNAFISH_props[0].intensity_mean
        if (total_int_DNAFISH > 1E8) & (total_int_DNAFISH < 1.1E8):
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
                    img_mCherry_pos_seg = ima.image_paste_to(img_mCherry_pos_seg, local_nuclear_seg_convex,
                                                                 [int(paste_to_centroid[0] - local_centroid[0]),
                                                                  int(paste_to_centroid[1] - local_centroid[1])])
                    img_mCherry_pos_ecseg = ima.image_paste_to(img_mCherry_pos_ecseg, local_ecseg,
                                                                 [int(paste_to_centroid[0] - local_centroid[0]),
                                                                  int(paste_to_centroid[1] - local_centroid[1])])
                    img_mCherry_pos_mCherry = ima.image_paste_to(img_mCherry_pos_mCherry, local_mCherry,
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
                    img_mCherry_neg_seg = ima.image_paste_to(img_mCherry_neg_seg, local_nuclear_seg_convex,
                                                             [int(paste_to_centroid[0] - local_centroid[0]),
                                                              int(paste_to_centroid[1] - local_centroid[1])])
                    img_mCherry_neg_ecseg = ima.image_paste_to(img_mCherry_neg_ecseg, local_ecseg,
                                                               [int(paste_to_centroid[0] - local_centroid[0]),
                                                                int(paste_to_centroid[1] - local_centroid[1])])
                    img_mCherry_neg_mCherry = ima.image_paste_to(img_mCherry_neg_mCherry, local_mCherry,
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
tif.imwrite("%s%s/grid/%s_pos_seg_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_seg)
tif.imwrite("%s%s/grid/%s_pos_ecseg_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_ecseg)
tif.imwrite("%s%s/grid/%s_neg_hoechst_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_hoechst)
tif.imwrite("%s%s/grid/%s_neg_DNAFISH_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_DNAFISH)
tif.imwrite("%s%s/grid/%s_neg_seg_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_seg)
tif.imwrite("%s%s/grid/%s_neg_ecseg_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_ecseg)
tif.imwrite("%s%s/grid/%s_pos_mCherry_%s.tif" % (output_dir, sample, sample, img_pos_count), img_mCherry_pos_mCherry)
tif.imwrite("%s%s/grid/%s_neg_mCherry_%s.tif" % (output_dir, sample, sample, img_neg_count), img_mCherry_neg_mCherry)

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_pos_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_pos_withmCh_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_pos_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
viewer.add_image(img_mCherry_pos_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
plt.imsave('%s%s/grid/%s_pos_seg_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_pos_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_pos_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_pos_hoechst_%s.tiff' % (output_dir, sample, sample, img_pos_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_neg_mCherry, blending='additive', colormap='red', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_neg_withmCh_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_neg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_seg, blending='additive', colormap='blue', contrast_limits=[0, 1])
viewer.add_image(img_mCherry_neg_ecseg, blending='additive', colormap='green', contrast_limits=[0, 1])
plt.imsave('%s%s/grid/%s_neg_seg_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_DNAFISH, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_neg_DNAFISH_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_mCherry_neg_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
plt.imsave('%s%s/grid/%s_neg_hoechst_%s.tiff' % (output_dir, sample, sample, img_neg_count), dis.blending(viewer))
viewer.close()