import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import shared.dataframe as dat
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif

# USAGE
# 01. if couldn't find automatically, need manual help
# 02. check about alignment

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8-1'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_fov_pixel = 58.6
refine_time = 3
refine_on_refine = 'YES'
# 766nm for Keyence 10x
# 360nm for 512x512 at 63x
# 58.6nm for 3144x3144 at 63x (0.0765)
# 180nm for 1024x1024 at 63x (0.2349)
# 720nm for 256x256 at 63x (0.9399)

# NO NEED TO CHANGE
print(sample)
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%salign/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
# treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s' % (sample, treatment, sample)
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s' % sample

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
total_fov = 16

img_before_GFP = skio.imread("%s/%s/%s_beforeFISH_GFP_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/%s/%s_beforeFISH_mCherry_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_before = img_before_GFP.copy()
img_before[img_before_mCherry > img_before_GFP] = img_before_mCherry[img_before_mCherry > img_before_GFP]
img_hoechst_local = np.zeros_like(img_after_hoechst)

if refine_on_refine == 'YES':
    topleft_final = pd.read_csv('%s/%s/align_topleft_final_refine_%s.txt' % (align_dir, sample, sample),
                                na_values=['.'], sep='\t')
else:
    topleft_final = pd.read_csv('%s/%s/align_topleft_final_manual_%s.txt' % (align_dir, sample, sample), na_values=['.'], sep='\t')
topleft_final['topleft_target'] = [dat.str_to_float(x) for x in topleft_final['topleft_target']]
topleft_final_refine = pd.DataFrame(columns=['sample', 'fov', 'topleft_target', 'min_ratio'])

for fov in range(total_fov):
    print(fov+1)
    img_hoechst = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov + 1), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize1 = cv2.flip(imutils.rotate(img_hoechst_resize, angle=-90), 0)
    topleft_target = topleft_final[(topleft_final['fov'] == fov+1)]['topleft_target'].tolist()[0]
    # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 5)
    for i in range(refine_time):
        topleft_target, min_ratio = ima.img_search_local(img_before, img_hoechst_resize1, topleft_target, 1)
    topleft_final_refine.loc[len(topleft_final_refine.index)] = [sample, fov + 1, topleft_target, min_ratio]
    _, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize1, [0, 0], topleft_target)
    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize_final, [0, 0])

    """viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()"""

topleft_final_refine.to_csv('%s/%s/align_topleft_final_refine_%s.txt' % (align_dir, sample, sample), index=False, sep='\t')
tif.imwrite("%s/%s/%s_fovFISH.tif" % (output_dir, sample, sample), img_hoechst_local)

viewer = napari.Viewer()
viewer.add_image(img_before, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()
# plt.imsave("%s/%s/%s_check_fovFISH-mergeFISH-afterFISH.tiff" % (output_dir, sample, sample), dis.blending(viewer))

print(sample)
print("step09 topleft individual FOV local search Done!")