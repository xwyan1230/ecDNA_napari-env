import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shared.image as ima
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_3_49pos'
dshape_factor = 0.415  # 58.7nm/142nm 0.413

topleft0 = [25300, 955]
n_col = 7

img_hoechst_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")
img_DNAFISH_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch00.tif" % (data_dir, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])

img_hoechst_IF = skio.imread("%sDM_324pos_merge/DM_hoechst_IF_22000-33601.tif" % data_dir, plugin="tifffile")
img_hoechst_local = np.zeros_like(img_hoechst_IF)

for fov in range(img_DNAFISH_stack.shape[0]):
    print(fov)
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_DNAFISH = img_DNAFISH_stack[fov, :, :]
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    s = img_hoechst_resize.shape[0]
    if fov == 0:
        topleft_target = topleft0
        data.loc[len(data.index)] = [sample, 0, topleft0[0], topleft0[1], 0]
    else:
        if fov//n_col == 0:
            topleft_target = [int(topleft0[0]-((fov%n_col)*1.1*s)), int(topleft0[1]+((fov//n_col)*1.1*s))]
        else:
            ref_fov = fov-(2*(fov%n_col)+1)
            topleft_target = [data[data['fov'] == ref_fov]['topleft_x'].tolist()[0], int(data[data['fov'] == ref_fov]['topleft_y'].tolist()[0]+1.1*s)]
            # topleft_target = [int(topleft0[0]-((n_col-1-fov%n_col)*1.1*s)), int(topleft0[1]+((fov//n_col)*1.1*s))]
        print(topleft_target)
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
                # print('%s/400: %s' % (i, min_ratio_temp))
                if min_ratio_temp < min_ratio:
                    min_ratio = min_ratio_temp
                    topleft_target = [topleft_ori1[0] + x * 5, topleft_ori1[1] + y * 5]

        data.loc[len(data.index)] = [sample, fov, topleft_target[0], topleft_target[1], min_ratio]
        print(min_ratio)
        print(topleft_target)

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize, [topleft_target[1], topleft_target[0]])

    img_hoechst_IF_cut = img_hoechst_IF.copy()
    img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]
    data.to_csv('%s%s_alignment.txt' % (output_dir, sample), index=False, sep='\t')

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_IF_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize, blending='additive', colormap='green', contrast_limits=[0, 45535])
    # viewer.add_image(minus, blending='additive')
    plt.imsave("%sDM_alignment/DM_DNAFISH_fov%s_alignment_local.tiff" % (output_dir, fov), dis.blending(viewer))
    viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_hoechst_IF, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 45535])
plt.imsave("%sDM_alignment/DM_DNAFISH_alignment_global.tiff" % output_dir, dis.blending(viewer))
napari.run()