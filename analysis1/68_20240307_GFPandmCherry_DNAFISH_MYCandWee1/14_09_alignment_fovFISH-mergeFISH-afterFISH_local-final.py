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

img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_hoechst_DNAFISH_cut = skio.imread("%s/%s/%s_hoechst_DNAFISH_withborder.tif" % (data_dir1, sample, sample), plugin='tifffile')
img_hoechst_local = np.zeros_like(img_after_hoechst)
img_screen_cut = skio.imread("%s/%s/%s_screen_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")

topleft_search = pd.read_csv('%s/%s/align_topleft_search_%s.txt' % (align_dir, sample, sample), na_values=['.'], sep='\t')
topleft_search['topleft_target'] = [dat.str_to_float(x) for x in topleft_search['topleft_target']]
topleft_final = pd.DataFrame(columns=['sample', 'fov', 'topleft_target_no', 'topleft_target', 'min_ratio'])
topleft_final_refine = pd.DataFrame(columns=['sample', 'fov', 'topleft_target', 'min_ratio'])

for i in range(total_fov):
    min_ratio_final = np.min(topleft_search[topleft_search['topleft_target_no'] == i]['min_ratio'].tolist())
    fov_final = topleft_search[(topleft_search['topleft_target_no'] == i) & (topleft_search['min_ratio'] == min_ratio_final)]['fov'].tolist()[0]
    topleft_target_final = topleft_search[(topleft_search['topleft_target_no'] == i) & (topleft_search['min_ratio'] == min_ratio_final)]['topleft_target'].tolist()[0]
    topleft_final.loc[len(topleft_final.index)] = [sample, fov_final, i, topleft_target_final, min_ratio_final]
topleft_final.to_csv('%s/%s/align_topleft_final_%s.txt' % (align_dir, sample, sample), index=False, sep='\t')

for fov in range(total_fov):
    print(fov+1)
    img_hoechst = skio.imread("%s/DNAFISH/%s/%s_%s_ch01.tif" % (data_dir, sample, filename, fov + 1), plugin="tifffile")
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst,
                                    dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize1 = cv2.flip(imutils.rotate(img_hoechst_resize, angle=-90), 0)
    if len(topleft_final[topleft_final['fov'] == fov+1]) != 0:
        min_ratio = np.min(topleft_final[topleft_final['fov'] == fov+1]['min_ratio'].tolist())
        topleft_target = topleft_final[(topleft_final['fov'] == fov+1) & (topleft_final['min_ratio'] == min_ratio)]['topleft_target'].tolist()[0]
        # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 20)
        # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 5)
        # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 1)
        topleft_final_refine.loc[len(topleft_final_refine.index)] = [sample, fov + 1, topleft_target, min_ratio]
        _, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize1, [0, 0], topleft_target)
        img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize_final, [0, 0])
    else:
        viewer = napari.Viewer()
        viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
        viewer.add_image(img_hoechst_resize1, blending='additive', colormap='green', contrast_limits=[0, 65535])
        viewer.add_image(img_screen_cut, blending='additive', contrast_limits=[0, img_screen_cut.max()])
        shapes = viewer.add_shapes(name='Shapes', ndim=2)
        napari.run()
        poly_data = shapes.data[0]

        topleft_target, min_ratio = ima.img_search_global(img_hoechst_DNAFISH_cut, img_hoechst_resize1, poly_data, 20)
        # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 5)
        # topleft_target, min_ratio = ima.img_search_local(img_after_hoechst, img_hoechst_resize1, topleft_target, 1)
        topleft_final_refine.loc[len(topleft_final_refine.index)] = [sample, fov + 1, topleft_target, min_ratio]
        _, img_hoechst_resize_final = ima.img_align_move(img_after_hoechst, img_hoechst_resize1, [0, 0], topleft_target)
        img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize_final, [0, 0])

    """viewer = napari.Viewer()
    viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_resize_final, blending='additive', colormap='green', contrast_limits=[0, 65535])
    viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()"""

topleft_final_refine.to_csv('%s/%s/align_topleft_final_manual_%s.txt' % (align_dir, sample, sample), index=False, sep='\t')
tif.imwrite("%s/%s/%s_fovFISH.tif" % (output_dir, sample, sample), img_hoechst_local)

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_DNAFISH_cut, blending='additive', contrast_limits=[0, 65535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 65535])
napari.run()
# plt.imsave("%s/%s/%s_check_fovFISH-mergeFISH-afterFISH.tiff" % (output_dir, sample, sample), dis.blending(viewer))

print(sample)
print("step09 topleft individual FOV local search Done!")