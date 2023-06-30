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

img_hoechst_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch00.tif" % (data_dir, sample, sample), plugin="tifffile")

for fov in range(img_DNAFISH_stack.shape[0]):
    print(fov)
    img_hoechst_IF = skio.imread("%sDM_324pos_merge/DM_hoechst_IF_22000-33601.tif" % data_dir, plugin="tifffile")
    print(img_hoechst_IF.shape)

    img_hoechst = img_hoechst_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
    print(img_hoechst_resize.shape)
    s = img_hoechst_resize.shape[0]
    print(s)

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    shapes_layer = viewer.layers['Shapes']
    poly_data = shapes.data[0]
    topleft = [poly_data[0][1], poly_data[0][0]]
    xrange = poly_data[1][1] - poly_data[0][1]
    yrange = poly_data[2][0] - poly_data[0][0]

    topleft = [25295, 950]
    """for x in range(xrange/50):
        for y in range(yrange/50):
            topleft_test = [topleft[0] + x*50, topleft[1] + y*50]"""

    img_hoechst_resize1 = np.concatenate([np.zeros(shape=[img_hoechst_resize.shape[0], topleft[0]]), img_hoechst_resize], axis=1)
    img_hoechst_resize2 = np.concatenate([np.zeros(shape=[topleft[1], img_hoechst_resize1.shape[1]]), img_hoechst_resize1], axis=0)

    img_hoechst_IF_cut = img_hoechst_IF.copy()
    img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft[1]):int(topleft[1]+s), int(topleft[0]):int(topleft[0]+s)]

    minus = img_hoechst_resize.astype(float) - img_hoechst_IF_cut.astype(float)
    minus[minus < 0] = 0
    min_ratio = minus.sum()/img_hoechst_resize.sum()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize, blending='additive', colormap='green', contrast_limits=[0, 45535])
    viewer.add_image(minus, blending='additive')
    # plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
    napari.run()

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize2, blending='additive', colormap='green', contrast_limits=[0, 45535])
    # plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
    napari.run()