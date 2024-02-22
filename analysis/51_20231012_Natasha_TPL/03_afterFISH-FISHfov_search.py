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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230724_analysis_Natasha_colcemid/TPL/mCh-DMSO_GFP-TPL/"
output_dir = "%s05_figures/" % master_folder

sample = 'mCh-DMSO_GFP-TPL'
interval = 50
dshape_factor = 0.471  # 766nm for Keyence 10x, 360nm for 512x512 at 63x, 58.6 for 3144x3144 at 63x (0.0765), 180 for 1024x1024 at 63x (0.2349)
total_fov = 49
n_col = 7
border = 0

data = pd.DataFrame(columns=['sample', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])

img_hoechst_DNAFISH_cut = skio.imread("%s/02_afterFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_overview_4_Merged_ch00.tif" % master_folder, plugin="tifffile")
ori_shape = img_hoechst_DNAFISH_cut.shape
img_hoechst_DNAFISH_cut = cv2.resize(img_hoechst_DNAFISH_cut, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                interpolation=cv2.INTER_AREA)
img_after_hoechst = skio.imread("%s/05_figures/01_alignment_global/4_img_hoechst_cut.tif" % master_folder, plugin="tifffile")
img_hoechst_local = np.zeros_like(img_after_hoechst)

fov = 0
filename = '230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1_s00'
img_hoechst = skio.imread("%s/03_DNAFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1/%s_ch00.tif" % (master_folder, filename), plugin="tifffile")
dshape_factor = 0.12
ori_shape = img_hoechst.shape
img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                interpolation=cv2.INTER_AREA)
s = img_hoechst_resize.shape[0]
print(s)
topleft_target = [img_after_hoechst.shape[0]-s-200-border, 250+border]
print(topleft_target)
topleft0, min_ratio0 = ima.img_search_local(img_after_hoechst, img_hoechst_resize, topleft_target, 5)
_, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize, [0, 0], topleft0)

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

for f in range(total_fov):
    fov = f
    print(fov)
    if fov < 10:
        filename = '230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1_s0%s' % fov
    else:
        filename = '230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1_s%s' % fov
    img_hoechst = skio.imread("%s/03_DNAFISH/230606_mChDMSO_GFPTPL_mChDMSO_GFPTPL_4_frame1/%s_ch00.tif" % (master_folder, filename), plugin="tifffile")
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
        topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize, topleft_target, 3)
        data.loc[len(data.index)] = [sample, fov, topleft_target[0], topleft_target[1], min_ratio]
        print(min_ratio)
        print(topleft_target)

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize, [int(topleft_target[1]), int(topleft_target[0])])
    if not os.path.exists("%s%s/" % (output_dir, sample)):
        os.makedirs("%s%s/" % (output_dir, sample))
    data.to_csv('%s%s/%s_4_alignment.txt' % (output_dir, sample, sample), index=False, sep='\t')

    """viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    # viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()"""

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
# viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave("%s/%s/%s_4_alignment_local.tiff" % (output_dir, sample, sample), dis.blending(viewer))
viewer.close()