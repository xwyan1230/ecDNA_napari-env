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
import nd2
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'F10'
total_fov = 16
dshape_factor = 0.145
pixel_size = 300/2720  # uM
align = pd.read_csv('%s/alignment.txt' % output_dir, na_values=['.'], sep='\t')
pd_loc = pd.read_csv("%s/%s/04_%s_location.txt" % (output_dir, sample, sample), na_values=['.'], sep='\t')
topleft_target = [align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_local')]['topleft_0'].tolist()[0],
                  align[(align['sample'] == sample) & (align['step'] == 'afterFISH-FISH_local')]['topleft_1'].tolist()[0]]
img_stack = nd2.imread('%sDNAFISH/%s.nd2' % (data_dir, sample))

data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y'])

img_before_GFP = skio.imread("%s/%s/03_%s_GFP_final.tif" % (output_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/%s/03_%s_mCherry_final.tif" % (output_dir, sample, sample), plugin="tifffile")
img_before = img_before_GFP.copy()
img_before[img_before_mCherry > img_before_GFP] = img_before_mCherry[img_before_mCherry > img_before_GFP]
img = imutils.rotate(img_before, angle=180)
"""img_after_hoechst = skio.imread("%s/%s/03_%s_hoechst_final.tif" % (output_dir, sample, sample), plugin="tifffile")
img = imutils.rotate(img_after_hoechst, angle=180)"""
img_hoechst_paste = np.zeros_like(img)

for fov in range(total_fov):
    img_hoechst = img_stack[fov, :, 0, :, :]
    img_hoechst_merge = img_hoechst.max(axis=0)
    ori_shape = img_hoechst_merge.shape
    img_hoechst_resize = cv2.resize(img_hoechst_merge, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    delta_x = int(((pd_loc['x'][fov] - min(pd_loc['x']))/pixel_size)*dshape_factor) + int(topleft_target[0])
    delta_y = -int(((pd_loc['y'][fov] - max(pd_loc['y']))/pixel_size)*dshape_factor) + int(topleft_target[1])

    _, img_hoechst_resize_final = ima.img_align_move(img, img_hoechst_resize, [0, 0], [delta_x, delta_y])

    viewer = napari.Viewer()
    viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 40000])
    # viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 20000])
    shapes = viewer.add_shapes(name='3nd fov check', ndim=2)
    napari.run()

    topleft0, min_ratio0 = ima.img_search_local(img, img_hoechst_resize, [delta_x, delta_y], 5)
    _, img_hoechst_resize_final = ima.img_align_move(img, img_hoechst_resize, [0, 0], topleft0)

    viewer = napari.Viewer()
    viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 40000])
    # viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 20000])
    shapes = viewer.add_shapes(name='3nd fov check', ndim=2)
    # napari.run()
    viewer.close()

    if len(shapes.data) == 0:
        data.loc[len(data.index)] = [sample, fov, topleft0[0], topleft0[1]]
        img_hoechst_paste = ima.image_paste_to(img_hoechst_paste, img_hoechst_resize,
                                               [int(topleft0[1]), int(topleft0[0])])
    """viewer = napari.Viewer()
    viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 40000])
    viewer.add_image(img_hoechst_paste, blending='additive', colormap='green', contrast_limits=[0, 20000])
    napari.run()"""

viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 40000])
viewer.add_image(img_hoechst_paste, blending='additive', colormap='green', contrast_limits=[0, 20000])
plt.imsave("%s/%s/08_%s_beforeFISH-FISHfov_global.tiff" % (output_dir, sample, sample), dis.blending(viewer))
napari.run()

data.to_csv('%s/%s/08_%s_alignment_global.txt' % (output_dir, sample, sample), index=False, sep='\t')

print("DONE!")