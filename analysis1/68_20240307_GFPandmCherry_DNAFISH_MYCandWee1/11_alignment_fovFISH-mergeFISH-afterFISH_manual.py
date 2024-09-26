import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
import os

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F4'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_fov_pixel = 180

# NO NEED TO CHANGE
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%salign/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s' % (sample, treatment, sample)

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
total_fov = 16
border = 50

data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])
img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_hoechst_DNAFISH_cut = skio.imread("%s/%s/%s_hoechst_DNAFISH_withborder.tif" % (data_dir1, sample, sample), plugin='tifffile')
img_hoechst_local = np.zeros_like(img_after_hoechst)
shape = img_after_hoechst.shape

img_screen = skio.imread("%s/screenshot/%s.png" % (data_dir, sample))
viewer = napari.Viewer()
viewer.add_image(img_screen, blending='additive', contrast_limits=[0, img_screen.max()])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]

img_screen_cut = img_screen[int(poly_data[0][0]):int(poly_data[2][0]), int(poly_data[0][1]):int(poly_data[1][1])]
ori_shape = img_screen_cut.shape
img_screen_cut_resize = cv2.resize(img_screen_cut, dsize=(shape[1]-2*border, shape[0]-2*border), interpolation=cv2.INTER_AREA)[:, :, 0]
_, img_screen_cut_resize = ima.img_align_move(img_after_hoechst, img_screen_cut_resize, [0, 0], [border, border])

for fov in range(total_fov):
    print(fov+1)
    img_hoechst = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov+1), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize1 = cv2.flip(imutils.rotate(img_hoechst_resize, angle=-90), 0)

    viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize1, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_screen_cut_resize, blending='additive', contrast_limits=[0, img_screen_cut_resize.max()])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
    poly_data = shapes.data[0]

    topleft_target, min_ratio = ima.img_search_global(img_hoechst_DNAFISH_cut, img_hoechst_resize1, poly_data, 20)
    topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 5)
    _, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize1, [0, 0], topleft_target)

    viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()

    data.loc[len(data.index)] = [sample, fov, topleft_target[0], topleft_target[1], min_ratio]

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize, [int(topleft_target[1]), int(topleft_target[0])])
    data.to_csv('%s%s/align_fovFISH-mergeFISH-afterFISH_%s.txt' % (output_dir, sample, sample), index=False, sep='\t')

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
# plt.imsave("%s/%s/%s_check_fovFISH-mergeFISH-afterFISH.tiff" % (output_dir, sample, sample), dis.blending(viewer))
viewer.close()