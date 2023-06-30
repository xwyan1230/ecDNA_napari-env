import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shared.image as ima
import tifffile as tif
from skimage.measure import label, regionprops
from skimage.draw import line_aa
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_3_49pos'
dshape_factor = 0.415  # 58.7nm/142nm 0.413

topleft0 = [25300, 955]
n_col = 7

img_hoechst_stack = skio.imread("%sDM_DNAFISH/%s/%s_RAW_ch01.tif" % (data_dir, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'FOV', 'index', 'centroid_nuclear'])
df_align = pd.read_csv('%s%s_alignment.txt' % (output_dir, sample), na_values=['.'], sep='\t')

img_hoechst_IF = skio.imread("%sDM_324pos_merge/DM_hoechst_IF_22000-33601.tif" % data_dir, plugin="tifffile")
img_hoechst_IF_seg = skio.imread("%sDM_324pos_merge/DM_nuclear_red_green_IF_22000-33601_seg_convex.tif" % data_dir, plugin="tifffile")

for fov in range(img_hoechst_stack.shape[0]):
    print("Analyzing %s, fov %s" % (sample, fov))
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_seg = skio.imread("%sDM_DNAFISH/%s/seg_tif/%s_%s_seg.tif" % (data_dir, sample, sample, fov), plugin="tifffile")
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(img_hoechst.shape[1] * dshape_factor), int(img_hoechst.shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)

    s = img_hoechst_resize.shape[0]
    img_alignment = np.zeros_like(img_hoechst_resize)
    topleft_target = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]

    img_hoechst_IF_cut = img_hoechst_IF.copy()
    img_hoechst_IF_cut = img_hoechst_IF_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]
    img_hoechst_IF_seg_cut = img_hoechst_IF_seg.copy()
    img_hoechst_IF_seg_cut = img_hoechst_IF_seg_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]

    DNAFISH_props = regionprops(img_seg)
    IF_props = regionprops(img_hoechst_IF_seg_cut)
    IF_centroid = [IF_props[x].centroid for x in range(len(IF_props))]

    for i in range(len(DNAFISH_props)):
        dis_square_min = 1000000
        neighbour_centriod = [0, 0]
        DNAFISH_centroid = DNAFISH_props[i].centroid
        DNAFISH_centroid_resize = [DNAFISH_centroid[0] * dshape_factor, DNAFISH_centroid[1] * dshape_factor]
        for j in range(len(IF_props)):
            dis_square_temp = (DNAFISH_centroid_resize[0] - IF_centroid[j][0])**2 + (DNAFISH_centroid_resize[1] - IF_centroid[j][1])**2
            if dis_square_temp < dis_square_min:
                dis_square_min = dis_square_temp
                neighbour_centriod = IF_centroid[j]
        if dis_square_min < 2500:
            neighbour_index = img_hoechst_IF_seg_cut[int(neighbour_centriod[0])][int(neighbour_centriod[1])]
            data.loc[len(data.index)] = [i, fov, neighbour_index, DNAFISH_centroid]
        else:
            data.loc[len(data.index)] = [i, fov, 0, DNAFISH_centroid]
            neighbour_centriod = DNAFISH_centroid_resize

        rr, cc, val = line_aa(int(DNAFISH_centroid_resize[0]), int(DNAFISH_centroid_resize[1]), int(neighbour_centriod[0]), int(neighbour_centriod[1]))
        img_alignment[rr, cc] = val * 65535

    tif.imwrite("%s/DM_DNAFISH/%s/seg_tif/%s_%s_alignment.tif" % (output_dir, sample, sample, fov), img_alignment)

    viewer = napari.Viewer()
    viewer.add_image(img_alignment, blending='additive', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_IF_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize, blending='additive', colormap='green', contrast_limits=[0, 45535])
    plt.imsave('%s/DM_DNAFISH/%s/color_img/%s_%s_img_alignment.tiff' % (output_dir, sample, sample, fov), dis.blending(viewer))
    viewer.close()

data.to_csv('%s%s_alignment_label.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")