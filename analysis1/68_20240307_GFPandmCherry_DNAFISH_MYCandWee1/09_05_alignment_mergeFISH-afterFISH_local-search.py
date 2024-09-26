import skimage.io as skio
import napari
import cv2
import shared.image as ima
import os
import shared.dataframe as dat
import pandas as pd
import tifffile as tif
import shared.display as dis
import matplotlib.pyplot as plt

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
output_dir = "%salign/" % master_folder
name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s_Merging_001' % sample
# treatment = name[name['sample'] == sample]['treatment'].tolist()[0]
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s_hoechst_Merging_001' % (sample, treatment,sample)
DNAFISH_hoechst_merge = "%s/DNAFISH/%s/%s_ch00.tif" % (data_dir, sample, filename)

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_merge_pixel * 1.0 / afterFISH_hoechst_pixel
interval = 5
align_global = pd.read_csv('%s/align_mergeFISH-afterFISH_global.txt' % align_dir, na_values=['.'], sep='\t')
align_global['global_topleft_target'] = [dat.str_to_float(x) for x in align_global['global_topleft_target']]
topleft_target = align_global[align_global['sample'] == sample]['global_topleft_target'].tolist()[0]
if os.path.exists("%s/align_mergeFISH-afterFISH_local.txt" % align_dir):
    align_local = pd.read_csv('%s/align_mergeFISH-afterFISH_local.txt' % align_dir, na_values=['.'], sep='\t')
else:
    align_local = pd.DataFrame(columns=['sample', 'local_topleft_target', 'local_min_ratio'])

# img
img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_final.tif" % (data_dir1, sample, sample), plugin="tifffile")
img = img_after_hoechst

# img_search
img_hoechst_DNAFISH = skio.imread(DNAFISH_hoechst_merge, plugin="tifffile")
ori_shape = img_hoechst_DNAFISH.shape
img_hoechst_DNAFISH1 = cv2.resize(img_hoechst_DNAFISH, dsize=(int(ori_shape[1]*dshape_factor),
                                                              int(ori_shape[0]*dshape_factor)),
                                  interpolation=cv2.INTER_AREA)
img_search = img_hoechst_DNAFISH1

# local_search
print("searching...")
# topleft_target, min_ratio = ima.img_search_local(img, img_search, topleft_target, 50)
topleft_target, min_ratio = ima.img_search_local(img, img_search, topleft_target, interval)

# local view
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]),
          int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave("%s%s/%s_check_mergeFISH-afterFISH_local.tiff" % (align_dir, sample, sample), dis.blending(viewer))
napari.run()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave("%s%s/%s_check_mergeFISH-afterFISH.tiff" % (align_dir, sample, sample), dis.blending(viewer))
napari.run()

tif.imwrite("%s/%s/%s_mergeFISH_final.tif" % (output_dir, sample, sample), img_search_final)

if sample in align_local['sample'].tolist():
    sample_index = align_local[align_local['sample'] == sample].index[0]
    align_local.loc[sample_index] = [sample, topleft_target, min_ratio]
else:
    align_local.loc[len(align_local.index)] = [sample, topleft_target, min_ratio]
align_local.to_csv('%s/align_mergeFISH-afterFISH_local.txt' % align_dir, index=False, sep='\t')

print(sample)
print(topleft_target)
print("step05 mergeFISH-afterFISH local search DONE!")