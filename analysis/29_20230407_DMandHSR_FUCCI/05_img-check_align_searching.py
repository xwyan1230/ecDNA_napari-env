import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_3_49pos'
dshape_factor = 0.415  # 58.7nm/142nm 0.413

topleft = [25295, 950]
# 25300, 955

img_hoechst_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch00.tif" % (data_dir, sample, sample), plugin="tifffile")

for fov in range(img_DNAFISH_stack.shape[0]):
    print(fov)
    img_hoechst_IF = skio.imread("%sDM_324pos_merge/DM_hoechst_IF_22000-33601.tif" % data_dir, plugin="tifffile")
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
    s = img_hoechst_resize.shape[0]

    """viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']
    poly_data = shapes.data[0]
    topleft_ori = [poly_data[0][1], poly_data[0][0]]
    xrange = poly_data[1][1] - poly_data[0][1]
    yrange = poly_data[2][0] - poly_data[0][0]
    print(xrange)
    print(yrange)

    # topleft = [25295, 950]
    min_ratio = 1
    topleft_target = topleft_ori
    i = 0
    for x in range(int(xrange/50)):
        for y in range(int(yrange/50)):
            i = i+1
            topleft = [topleft_ori[0] + x*50, topleft_ori[1] + y*50]
            img_hoechst_IF_cut = img_hoechst_IF.copy()
            img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft[1]):int(topleft[1]+s), int(topleft[0]):int(topleft[0]+s)]

            minus = img_hoechst_resize.astype(float) - img_hoechst_IF_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum()/img_hoechst_resize.sum()
            print('%s/%s: %s' % (i, int(xrange/50) * int(yrange/50), min_ratio_temp))
            if min_ratio_temp < min_ratio:
                min_ratio = min_ratio_temp
                topleft_target = [topleft_ori[0] + x*50, topleft_ori[1] + y*50]
    print(min_ratio)
    print(topleft_target)"""

    topleft_target = [25315, 945]
    min_ratio = 1
    topleft_ori1 = [topleft_target[0] - 50, topleft_target[1] - 50]
    i = 0
    for x in range(20):
        for y in range(20):
            i = i+1
            topleft = [topleft_ori1[0] + x * 5, topleft_ori1[1] + y * 5]
            img_hoechst_IF_cut = img_hoechst_IF.copy()
            img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft[1]):int(topleft[1] + s),
                                 int(topleft[0]):int(topleft[0] + s)]

            minus = img_hoechst_resize.astype(float) - img_hoechst_IF_cut.astype(float)
            minus[minus < 0] = 0
            min_ratio_temp = minus.sum() / img_hoechst_resize.sum()
            print('%s/400: %s' % (i, min_ratio_temp))
            if min_ratio_temp < min_ratio:
                min_ratio = min_ratio_temp
                topleft_target = [topleft_ori1[0] + x * 5, topleft_ori1[1] + y * 5]

    print(min_ratio)
    print(topleft_target)

    img_hoechst_resize1 = np.concatenate(
        [np.zeros(shape=[img_hoechst_resize.shape[0], int(topleft_target[0])]), img_hoechst_resize], axis=1)
    img_hoechst_resize2 = np.concatenate(
        [np.zeros(shape=[int(topleft_target[1]), img_hoechst_resize1.shape[1]]), img_hoechst_resize1], axis=0)

    img_hoechst_IF_cut = img_hoechst_IF.copy()
    img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize, blending='additive', colormap='green', contrast_limits=[0, 45535])
    # viewer.add_image(minus, blending='additive')
    plt.imsave("%sDM_alignment/DM_DNAFISH_fov0_alignment_local.tiff" % output_dir, dis.blending(viewer))
    viewer.close()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize2, blending='additive', colormap='green', contrast_limits=[0, 45535])
    plt.imsave("%sDM_alignment/DM_DNAFISH_fov0_alignment_global.tiff" % output_dir, dis.blending(viewer))
    napari.run()