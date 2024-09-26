import skimage.io as skio
import napari
import imutils
import pandas as pd
import os
import shared.image as ima
import numpy as np

# USAGE
# 01. square the region that you think where the topleft corner of the given image (red, part in green image) in the
# original image (blue)
# 02. check about alignment

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F5'
search_offset = 2000

# NO NEED TO CHANGE
print(sample)
data_dir = "%sdata/" % master_folder
output_dir = "%salign/" % master_folder
beforeFISH_GFP = "%s/beforeFISH/%s_GFP.tif" % (data_dir, sample)
beforeFISH_mCherry = "%s/beforeFISH/%s_mCherry.tif" % (data_dir, sample)
afterFISH_hoechst = "%s/afterFISH/%s_hoechst.tif" % (data_dir, sample)

# DO NOT CHANGE
interval = 50
if os.path.exists("%s/align_beforeFISH-afterFISH_global.txt" % output_dir):
    align_global = pd.read_csv('%s/align_beforeFISH-afterFISH_global.txt' % output_dir, na_values=['.'], sep='\t')
else:
    align_global = pd.DataFrame(columns=['sample', 'global_topleft_target', 'global_min_ratio', 'search_offset'])

# img
img_before_GFP = skio.imread(beforeFISH_GFP, plugin="tifffile")[:, :, 1]
img_before_mCherry = skio.imread(beforeFISH_mCherry, plugin="tifffile")[:, :, 0]
img_before = img_before_GFP.copy()
img_before[img_before_mCherry > img_before_GFP] = img_before_mCherry[img_before_mCherry > img_before_GFP]
img_before1 = imutils.rotate(img_before, angle=90)
img = img_before1

# img_search
img_after_hoechst = skio.imread(afterFISH_hoechst, plugin="tifffile")[:, :, 2]
img_after_hoechst1 = imutils.rotate(img_after_hoechst, angle=90)
img_after_hoechst1_test = img_after_hoechst1[search_offset:search_offset+3000, search_offset:search_offset+3000]
img_after_hoechst1_test1 = np.concatenate([np.zeros(shape=[search_offset, img_after_hoechst1_test.shape[1]]),
                                           img_after_hoechst1_test], axis=0)
img_after_hoechst1_test2 = np.concatenate([np.zeros(shape=[img_after_hoechst1_test1.shape[0], search_offset]),
                                           img_after_hoechst1_test1], axis=1)
img_search = img_after_hoechst1_test

# search
viewer = napari.Viewer()
viewer.add_image(img, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst1, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst1_test2, blending='additive', colormap='red', contrast_limits=[0, 65535])
shapes = viewer.add_shapes(name='Shapes', ndim=2)
napari.run()
poly_data = shapes.data[0]
topleft_target, min_ratio = ima.img_search_global(img, img_search, poly_data, interval)

# local view
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

if sample in align_global['sample'].tolist():
    sample_index = align_global[align_global['sample'] == sample].index[0]
    align_global.loc[sample_index] = [sample, topleft_target, min_ratio, search_offset]
else:
    align_global.loc[len(align_global.index)] = [sample, topleft_target, min_ratio, search_offset]
if not os.path.exists("%s/" % output_dir):
    os.makedirs("%s/" % output_dir)
align_global.to_csv('%s/align_beforeFISH-afterFISH_global.txt' % output_dir, index=False, sep='\t')

print(sample)
print(topleft_target)
print(search_offset)
print("step01 beforeFISH-afterFISH global search DONE!")
