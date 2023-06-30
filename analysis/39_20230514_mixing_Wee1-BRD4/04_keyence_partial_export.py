import skimage.io as skio
import napari
import cv2
import imutils
import shared.display as dis
import shared.image as ima
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
dshape_factor = 0.47  # 766nm for Keyence 10x, 360nm for 512x512 at 63x
topleft_target = [3865.6283976329314, 2319.720277882853]
border = 200

# topleft_target = [3313.627431155869, 1430.784821051334]  2_Ctrl-GFP-6hr_Ctrl-mCh-6hr [3333.627431155869, 1430.784821051334] 1.3
# topleft_target = [3753.8435086956088, 5124.393828360026]  5_Ctrl-GFP-6hr_Wee1-mCh-6hr [3733.8435086956088, 5094.393828360026] 1.1
# topleft_target = [5187.786656194921, 4282.431129459471]  6_Wee1-GFP-6hr_Ctrl-mCh-6hr [5177.786656194921, 4282.431129459471] 1.1
# topleft_target = [2615.5115503413313, 4631.470673742484]  9_Ctrl-GFP-6hr_BRD4-mCh-6hr [2635.5115503413313, 4641.470673742484] 0.5
# topleft_target = [3875.6283976329314, 2309.720277882853]  10_BRD4-GFP-6hr_Ctrl-mCh-6hr [3865.6283976329314, 2319.720277882853] 0.3

img_hoechst_Keyence = skio.imread("%s%s/Hoechst.tif" % (data_dir, sample), plugin="tifffile")[:, :, 2]
img_mCherry_Keyence = skio.imread("%s%s/mCherry.tif" % (data_dir, sample), plugin="tifffile")[:, :, 0]
img_GFP_Keyence = skio.imread("%s%s/GFP.tif" % (data_dir, sample), plugin="tifffile")[:, :, 1]
img_hoechst_DNAFISH = skio.imread("%s%s/%s_%s_hoechst_512_Merging_001_RAW_ch00.tif" % (data_dir, sample, exp, sample), plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor), int(ori_shape[0]*dshape_factor)), interpolation=cv2.INTER_AREA)
img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft_target[1])-border:int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0])+border, int(topleft_target[0])-border:int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])+border]
img_mCherry_Keyence_cut = img_mCherry_Keyence.copy()
img_mCherry_Keyence_cut = img_mCherry_Keyence_cut[int(topleft_target[1])-border:int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0])+border, int(topleft_target[0])-border:int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])+border]
img_GFP_Keyence_cut = img_GFP_Keyence.copy()
img_GFP_Keyence_cut = img_GFP_Keyence_cut[int(topleft_target[1])-border:int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0])+border, int(topleft_target[0])-border:int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])+border]
tif.imwrite("%s%s/Hoechst_cut.tif" % (output_dir, sample), img_hoechst_Keyence_cut)
tif.imwrite("%s%s/mCherry_cut.tif" % (output_dir, sample), img_mCherry_Keyence_cut)
tif.imwrite("%s%s/GFP_cut.tif" % (output_dir, sample), img_GFP_Keyence_cut)
print("DONE!")