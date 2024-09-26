import skimage.io as skio
import napari
import pandas as pd
import shared.dataframe as dat
import shared.image as ima
import os
import numpy as np
import tifffile as tif

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
align_local = pd.read_csv('%s/align_beforeFISH-afterFISH_local.txt' % align_dir, na_values=['.'], sep='\t')
align_local['local_topleft_target'] = [dat.str_to_float(x) for x in align_local['local_topleft_target']]
topleft_target = align_local[align_local['sample'] == sample]['local_topleft_target'].tolist()[0]
search_offset = align_local[align_local['sample'] == sample]['search_offset'].tolist()[0]
topleft_target = list(np.array(topleft_target) - np.array([search_offset, search_offset]))
print(topleft_target)  # positive change 2nd image, negative change 1st image

img_before_GFP = skio.imread(beforeFISH_GFP, plugin="tifffile")[:, :, 1]
img_before_mCherry = skio.imread(beforeFISH_mCherry, plugin="tifffile")[:, :, 0]
img_after_hoechst = skio.imread(afterFISH_hoechst, plugin="tifffile")[:, :, 2]
img_before_GFP_final, img_after_hoechst_final = ima.img_align_move(img_before_GFP, img_after_hoechst, [90, 90],
                                                                   topleft_target)
img_before_mCherry_final, _ = ima.img_align_move(img_before_mCherry, img_after_hoechst, [90, 90], topleft_target)

if not os.path.exists("%s%s/" % (output_dir, sample)):
    os.makedirs("%s%s/" % (output_dir, sample))
tif.imwrite("%s/%s/%s_beforeFISH_GFP_final.tif" % (output_dir, sample, sample), img_before_GFP_final)
tif.imwrite("%s/%s/%s_beforeFISH_mCherry_final.tif" % (output_dir, sample, sample), img_before_mCherry_final)
tif.imwrite("%s/%s/%s_afterFISH_hoechst_final.tif" % (output_dir, sample, sample), img_after_hoechst_final)

viewer = napari.Viewer()
viewer.add_image(img_before_GFP_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
viewer.add_image(img_before_mCherry_final, blending='additive', colormap='red', contrast_limits=[0, 65535])
viewer.add_image(img_after_hoechst_final, blending='additive', colormap='blue', contrast_limits=[0, 65535])
# plt.imsave("%s%s/DM_alignment.tiff" % (output_dir, sample), dis.blending(viewer))
napari.run()

print(sample)
print("step03 beforeFISH-afterFISH image check and save DONE!")