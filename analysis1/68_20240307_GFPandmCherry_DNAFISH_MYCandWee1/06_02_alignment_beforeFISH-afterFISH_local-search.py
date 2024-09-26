import skimage.io as skio
import napari
import imutils
import shared.dataframe as dat
import pandas as pd
import shared.image as ima
import shared.display as dis
import matplotlib.pyplot as plt
import os

# USAGE
# 01. no need to do anything
# 02. check about alignment

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F5'

# NO NEED TO CHANGE
print(sample)
data_dir = "%sdata/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%sprocessed/" % master_folder
beforeFISH_GFP = "%s/beforeFISH/%s_GFP.tif" % (data_dir, sample)
beforeFISH_mCherry = "%s/beforeFISH/%s_mCherry.tif" % (data_dir, sample)
afterFISH_hoechst = "%s/afterFISH/%s_hoechst.tif" % (data_dir, sample)

# DO NOT CHANGE
align_global = pd.read_csv('%s/align_beforeFISH-afterFISH_global.txt' % align_dir, na_values=['.'], sep='\t')
align_global['global_topleft_target'] = [dat.str_to_float(x) for x in align_global['global_topleft_target']]
topleft_target = align_global[align_global['sample'] == sample]['global_topleft_target'].tolist()[0]
search_offset = align_global[align_global['sample'] == sample]['search_offset'].tolist()[0]
interval = 5
if os.path.exists("%s/align_beforeFISH-afterFISH_local.txt" % align_dir):
    align_local = pd.read_csv('%s/align_beforeFISH-afterFISH_local.txt' % align_dir, na_values=['.'], sep='\t')
else:
    align_local = pd.DataFrame(columns=['sample', 'local_topleft_target', 'local_min_ratio', 'search_offset'])

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
img_search = img_after_hoechst1_test

# local_search
print("searching...")
topleft_target, min_ratio = ima.img_search_local(img, img_search, topleft_target, interval)

# local view
img_cut = img.copy()
img_cut = img_cut[int(topleft_target[1]):int(topleft_target[1]+img_search.shape[0]),
          int(topleft_target[0]):int(topleft_target[0]+img_search.shape[1])]

viewer = napari.Viewer()
viewer.add_image(img_cut, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search, blending='additive', colormap='green', contrast_limits=[0, 65535])
if not os.path.exists("%s%s/" % (align_dir, sample)):
    os.makedirs("%s%s/" % (align_dir, sample))
plt.imsave("%s%s/%s_check_before-afterFISH_local.tiff" % (align_dir, sample, sample), dis.blending(viewer))
napari.run()

# global view
img_final, img_search_final = ima.img_align_move(img, img_search, [0, 0], topleft_target)

viewer = napari.Viewer()
viewer.add_image(img_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_search_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
plt.imsave("%s%s/%s_check_before-afterFISH.tiff" % (align_dir, sample, sample), dis.blending(viewer))
napari.run()

if sample in align_local['sample'].tolist():
    sample_index = align_local[align_local['sample'] == sample].index[0]
    align_local.loc[sample_index] = [sample, topleft_target, min_ratio, search_offset]
else:
    align_local.loc[len(align_local.index)] = [sample, topleft_target, min_ratio, search_offset]
align_local.to_csv('%s/align_beforeFISH-afterFISH_local.txt' % align_dir, index=False, sep='\t')

print(sample)
print(topleft_target)
print(search_offset)
print("step02 beforeFISH-afterFISH local search DONE!")
