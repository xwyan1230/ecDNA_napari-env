import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.image as ima
import tifffile as tif
from skimage.measure import label, regionprops
from skimage.draw import line_aa
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'
batch = 2
dshape_factor = 0.0765
a_target = 0.3

n_col = 10

img_hoechst_stack = skio.imread("%s%s/%s_%s_%s_50pos_RAW_ch00.tif" % (data_dir, sample, exp, sample, batch), plugin="tifffile")

data = pd.DataFrame(columns=['nuclear', 'batch', 'FOV', 'index', 'centroid_nuclear'])
df_align = pd.read_csv('%s%s/%s_%s_alignment.txt' % (data_dir1, sample, sample, batch), na_values=['.'], sep='\t')

img_hoechst_Keyence = skio.imread("%s%s/Hoechst_cut.tif" % (data_dir, sample), plugin="tifffile")
img_hoechst_Keyence_seg = skio.imread("%s%s/seg_keyence/Hoechst_keyence_seg_convex.tif" % (data_dir1, sample), plugin="tifffile")

for fov in range(img_hoechst_stack.shape[0]):
# for fov in range(60):
    print("Analyzing %s, fov %s" % (sample, fov))
    img_hoechst = img_hoechst_stack[fov, :, :]
    img_seg = (skio.imread("%s%s/seg_tif/%s_%s_%s_seg.tif" % (data_dir1, sample, sample, batch, fov), plugin="tifffile")).astype(np.int16)
    img_seg_rotate = imutils.rotate(img_seg, angle=a_target)
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(img_hoechst.shape[1] * dshape_factor), int(img_hoechst.shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize_rotate = imutils.rotate(img_hoechst_resize, angle=a_target)
    s = img_hoechst_resize_rotate.shape[0]

    img_alignment = np.zeros_like(img_hoechst_resize_rotate)
    topleft_target = [df_align['topleft_x'][fov], df_align['topleft_y'][fov]]

    img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
    img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]
    img_hoechst_Keyence_seg_cut = img_hoechst_Keyence_seg.copy()
    img_hoechst_Keyence_seg_cut = img_hoechst_Keyence_seg_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]

    DNAFISH_props = regionprops(img_seg_rotate)
    Keyence_props = regionprops(img_hoechst_Keyence_seg_cut)
    Keyence_centroid = [Keyence_props[x].centroid for x in range(len(Keyence_props))]

    for i in range(len(DNAFISH_props)):
        dis_square_min = 1000000
        neighbour_centriod = [0, 0]
        DNAFISH_centroid = DNAFISH_props[i].centroid
        DNAFISH_centroid_resize = [DNAFISH_centroid[0] * dshape_factor, DNAFISH_centroid[1] * dshape_factor]
        for j in range(len(Keyence_props)):
            dis_square_temp = (DNAFISH_centroid_resize[0] - Keyence_centroid[j][0]) ** 2 + (DNAFISH_centroid_resize[1] - Keyence_centroid[j][1]) ** 2
            if dis_square_temp < dis_square_min:
                dis_square_min = dis_square_temp
                neighbour_centriod = Keyence_centroid[j]
        if dis_square_min < 500:
            neighbour_index = img_hoechst_Keyence_seg_cut[int(neighbour_centriod[0])][int(neighbour_centriod[1])]
            data.loc[len(data.index)] = [i, batch, fov, neighbour_index, DNAFISH_centroid]
        else:
            data.loc[len(data.index)] = [i, batch, fov, 0, DNAFISH_centroid]
            neighbour_centriod = DNAFISH_centroid_resize

        rr, cc, val = line_aa(int(DNAFISH_centroid_resize[0]), int(DNAFISH_centroid_resize[1]), int(neighbour_centriod[0]), int(neighbour_centriod[1]))
        img_alignment[rr, cc] = val * 65535
    if not os.path.exists("%s%s/align_assign/" % (output_dir, sample)):
        os.makedirs("%s%s/align_assign/" % (output_dir, sample))
    tif.imwrite("%s%s/align_assign/%s_%s_%s_alignment.tif" % (output_dir, sample, sample, batch, fov), img_alignment)

    viewer = napari.Viewer()
    viewer.add_image(img_alignment, blending='additive', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_Keyence_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize_rotate, blending='additive', colormap='green', contrast_limits=[0, 45535])
    plt.imsave('%s%s/align_assign/%s_%s_%s_img_alignment.tiff' % (output_dir, sample, sample, batch, fov), dis.blending(viewer))
    viewer.close()

data.to_csv('%s%s/alignment_label_%s.txt' % (output_dir, sample, batch), index=False, sep='\t')

print("DONE!")