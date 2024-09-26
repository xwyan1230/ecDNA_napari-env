import skimage.io as skio
import napari
import shared.image as ima
import cv2
import pandas as pd
import shared.dataframe as dat
import tifffile as tif

# USAGE
# 01. no need to do anything
# 02. check about alignment

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8-1'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_merge_pixel = 720
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
output_dir = "%sprocessed/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
# treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s_hoechst_Merging_001' % (sample, treatment, sample)
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s_Merging_001' % sample
DNAFISH_hoechst_merge = "%s/DNAFISH/%s/%s_ch00.tif" % (data_dir, sample, filename)

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_merge_pixel * 1.0 / afterFISH_hoechst_pixel
align_local = pd.read_csv('%s/align_mergeFISH-afterFISH_local.txt' % align_dir, na_values=['.'], sep='\t')
align_local['local_topleft_target'] = [dat.str_to_float(x) for x in align_local['local_topleft_target']]
topleft_target = align_local[align_local['sample'] == sample]['local_topleft_target'].tolist()[0]
border = 100

img_before_GFP = skio.imread("%s/%s/%s_beforeFISH_GFP_final.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_before_mCherry = skio.imread("%s/%s/%s_beforeFISH_mCherry_final.tif" % (data_dir1, sample, sample),
                                 plugin="tifffile")
img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_final.tif" % (data_dir1, sample, sample), plugin="tifffile")

img_hoechst_DNAFISH = skio.imread(DNAFISH_hoechst_merge, plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor),
                                                              int(ori_shape[0]*dshape_factor)),
                                  interpolation=cv2.INTER_AREA)

range1 = int(topleft_target[1])-border
range2 = int(topleft_target[1]+img_hoechst_DNAFISH1.shape[0])+border
range3 = int(topleft_target[0])-border
range4 = int(topleft_target[0]+img_hoechst_DNAFISH1.shape[1])+border
img_before_GFP_cut = img_before_GFP[range1:range2, range3:range4]
img_before_mCherry_cut = img_before_mCherry[range1:range2, range3:range4]
img_after_hoechst_cut = img_after_hoechst[range1:range2, range3:range4]
_, img_hoechst_DNAFISH_final = ima.img_align_move(img_after_hoechst, img_hoechst_DNAFISH1, [0, 0], [border, border])

tif.imwrite("%s/%s/%s_beforeFISH_GFP_cut.tif" % (output_dir, sample, sample), img_before_GFP_cut)
tif.imwrite("%s/%s/%s_beforeFISH_mCherry_cut.tif" % (output_dir, sample, sample), img_before_mCherry_cut)
tif.imwrite("%s/%s/%s_afterFISH_hoechst_cut.tif" % (output_dir, sample, sample), img_after_hoechst_cut)
tif.imwrite("%s/%s/%s_hoechst_DNAFISH_withborder.tif" % (output_dir, sample, sample), img_hoechst_DNAFISH_final)

viewer = napari.Viewer()
viewer.add_image(img_before_GFP_cut, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry_cut, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_final, blending='additive', colormap='blue', contrast_limits=[0, 45535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()

print(sample)
print("step06 mergeFISH-afterFISH image check and save DONE!")