import skimage.io as skio
import napari
import shared.image as ima
import cv2
import tifffile as tif
import imutils
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

row = 'G'
sample = 'G9'
dshape_factor = 0.2349  # 766nm for Keyence 10x, 360nm for 512x512 at 63x, 58.6 for 3144x3144 at 63x (0.0765), 180 for 1024x1024 at 63x (0.2349)
topleft_target = [3311.4535513688943, 6905.216405179622]
border = 50

img_before_GFP = skio.imread("%s/beforeFISH/%s/%s_GFP_final.tif" % (data_dir, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/beforeFISH/%s/%s_mCherry_final.tif" % (data_dir, sample, sample), plugin="tifffile")
img_after_hoechst = skio.imread("%s/afterFISH/%s/%s_hoechst_final.tif" % (data_dir, sample, sample), plugin="tifffile")

# filename = '20230601_CRISPRko_48hr_DNAFISH_%s_%s_scan_1024_Merging_001_ch00' % (row, sample)
filename = '20230602_DM_CRISPRko_48hr_DNAFISH_%s_%s_scan_1024_Merging_001_ch00' % (row, sample)
img_hoechst_DNAFISH = skio.imread("%s/FISH/%s/%s.tif" % (data_dir, sample, filename), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)

range1 = int(topleft_target[1])-border
range2 = int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0])+border
range3 = int(topleft_target[0])-border
range4 = int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])+border
img_before_GFP_cut = img_before_GFP[range1:range2, range3:range4]
img_before_mCherry_cut = img_before_mCherry[range1:range2, range3:range4]
img_after_hoechst_cut = img_after_hoechst[range1:range2, range3:range4]
"""ori_shape = img_before_GFP_cut.shape
img_before_GFP_cut_enlarge = cv2.resize(img_before_GFP_cut, dsize=(int(ori_shape[1]*(1/dshape_factor)), int(ori_shape[0]*(1/dshape_factor))), interpolation=cv2.INTER_AREA)
img_before_mCherry_cut_enlarge = cv2.resize(img_before_mCherry_cut, dsize=(int(ori_shape[1]*(1/dshape_factor)), int(ori_shape[0]*(1/dshape_factor))), interpolation=cv2.INTER_AREA)
img_after_hoechst_cut_enlarge = cv2.resize(img_after_hoechst_cut, dsize=(int(ori_shape[1]*(1/dshape_factor)), int(ori_shape[0]*(1/dshape_factor))), interpolation=cv2.INTER_AREA)"""
_, img_hoechst_DNAFISH_final = ima.img_align_move(img_after_hoechst, img_hoechst_DNAFISH1, [0, 0], [border, border])

tif.imwrite("%s/beforeFISH/%s/%s_GFP_cut.tif" % (output_dir, sample, sample), img_before_GFP_cut)
tif.imwrite("%s/beforeFISH/%s/%s_mCherry_cut.tif" % (output_dir, sample, sample), img_before_mCherry_cut)
tif.imwrite("%s/afterFISH/%s/%s_hoechst_cut.tif" % (output_dir, sample, sample), img_after_hoechst_cut)
tif.imwrite("%s/FISH/%s/%s_hoechst_DNAFISH_withborder.tif" % (output_dir, sample, sample), img_hoechst_DNAFISH_final)

viewer = napari.Viewer()
viewer.add_image(img_before_GFP_cut, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry_cut, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_final, blending='additive', colormap='blue', contrast_limits=[0, 45535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()