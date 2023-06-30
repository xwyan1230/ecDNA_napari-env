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
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%salignment/" % master_folder

row = 'C'
sample = 'C6'
batch = 2
total_fov = 12
dshape_factor = 0.0765
n_col = 4
border = 50

data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])

img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst_cut.tif" % (data_dir, sample, sample), plugin="tifffile")
img_hoechst_DNAFISH_cut = skio.imread("%s/FISH/%s/%s_hoechst_DNAFISH_withborder.tif" % (data_dir, sample, sample), plugin='tifffile')
img_hoechst_local = np.zeros_like(img_after_hoechst)

fov = 0
filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s00_ch01' % (row, sample)
img_hoechst = skio.imread("%s/FISH/%s/%s.tif" % (data_dir, sample, filename), plugin="tifffile")
ori_shape = img_hoechst.shape
img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                interpolation=cv2.INTER_AREA)
s = img_hoechst_resize.shape[0]
print(s)
topleft_target = [img_after_hoechst.shape[0]-s-250-border, 250+border+1.1*s]
print(topleft_target)
topleft0, min_ratio0 = ima.img_search_local(img_after_hoechst, img_hoechst_resize, topleft_target, 5)
_, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize, [0, 0], topleft0)

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

for fov in range(total_fov):
    print(fov)
    if fov < 10:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s0%s_ch01' % (row, sample, fov)
    else:
        filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_12pos_s%s_ch01' % (row, sample, fov)
    img_hoechst = skio.imread("%s/FISH/%s/%s.tif" % (data_dir, sample, filename), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    s = img_hoechst_resize.shape[0]
    if fov == 0:
        topleft_target = topleft0
        data.loc[len(data.index)] = [sample, 0, topleft0[0], topleft0[1], min_ratio0]
    else:
        if fov//n_col == 0:
            topleft_target = [int(topleft0[0]-((fov%n_col)*1.1*s)), int(topleft0[1]+((fov//n_col)*1.1*s))]
        else:
            ref_fov = fov-(2*(fov%n_col)+1)
            topleft_target = [data[data['fov'] == ref_fov]['topleft_x'].tolist()[0], int(data[data['fov'] == ref_fov]['topleft_y'].tolist()[0]+1.1*s)]
        print(topleft_target)
        topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize, topleft_target, 5)
        data.loc[len(data.index)] = [sample, fov, topleft_target[0], topleft_target[1], min_ratio]
        print(min_ratio)
        print(topleft_target)

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize, [int(topleft_target[1]), int(topleft_target[0])])
    data.to_csv('%s/%s/%s_alignment_%s.txt' % (output_dir, sample, sample, batch), index=False, sep='\t')

    """viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()"""

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
# viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()
plt.imsave("%s/%s/%s_alignment_local_%s.tiff" % (output_dir, sample, sample, batch), dis.blending(viewer))
viewer.close()