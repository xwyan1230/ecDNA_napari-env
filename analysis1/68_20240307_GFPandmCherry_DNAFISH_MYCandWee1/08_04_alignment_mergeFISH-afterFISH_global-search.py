import skimage.io as skio
import napari
import cv2
import pandas as pd
import shared.image as ima
import tifffile as tif
import os

# USAGE
# 01. square the region that you think where the topleft corner of the search image (white in screen image) in the
# original image (blue)
# 02. check about the alignment

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F8-1'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_merge_pixel = 720
ref_sample = 'F4'
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
# filename = '20240301_sp8_GFPandmCherry_MYCandWee1_4day_DNAFISH_%s_%s_%s_hoechst_Merging_001' % (sample, treatment,sample)
filename = '20240311_GFPandmCherry_DNAFISH_F8_re-capture_%s_Merging_001' % sample
DNAFISH_hoechst_merge = "%s/DNAFISH/%s/%s_ch00.tif" % (data_dir, sample, filename)

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_merge_pixel * 1.0 / afterFISH_hoechst_pixel
interval = 50
if os.path.exists("%s/align_mergeFISH-afterFISH_global.txt" % output_dir):
    align_global = pd.read_csv('%s/align_mergeFISH-afterFISH_global.txt' % output_dir, na_values=['.'], sep='\t')
else:
    align_global = pd.DataFrame(columns=['sample', 'global_topleft_target', 'global_min_ratio'])

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

# referene images
img_screen = skio.imread("%s/screenshot/%s.png" % (data_dir, sample))
img_screen1 = cv2.resize(img_screen, dsize=(9500, 9500), interpolation=cv2.INTER_AREA)

if os.path.exists("%s/%s/%s_mergeFISH_final.tif" % (align_dir, ref_sample, ref_sample)):
    img_ref = skio.imread("%s/%s/%s_mergeFISH_final.tif" % (align_dir, ref_sample, ref_sample))
    img_ref_screen = skio.imread("%s/screenshot/%s.png" % (data_dir, ref_sample))
    img_ref_screen1 = cv2.resize(img_ref_screen, dsize=(9500, 9500), interpolation=cv2.INTER_AREA)
    # global search
    viewer = napari.Viewer()
    viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    # viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_screen1, blending='additive', contrast_limits=[0, img_screen.max()])
    viewer.add_image(img_ref, blending='additive', contrast_limits=[0, img_ref.max()])
    viewer.add_image(img_ref_screen1, blending='additive', contrast_limits=[0, img_ref_screen.max()])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
else:
    viewer = napari.Viewer()
    viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_screen1, blending='additive', contrast_limits=[0, img_screen.max()])
    shapes = viewer.add_shapes(name='Shapes', ndim=2)
    napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)

# local view
print(topleft_target)
print(sample)
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]),
          int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()

tif.imwrite("%s/%s/%s_mergeFISH_final.tif" % (output_dir, sample, sample), img_search_final)

if sample in align_global['sample'].tolist():
    sample_index = align_global[align_global['sample'] == sample].index[0]
    align_global.loc[sample_index] = [sample, topleft_target, min_ratio]
else:
    align_global.loc[len(align_global.index)] = [sample, topleft_target, min_ratio]
align_global.to_csv('%s/align_mergeFISH-afterFISH_global.txt' % output_dir, index=False, sep='\t')

print(sample)
print(topleft_target)
print("step04 mergeFISH-afterFISH global search DONE!")

