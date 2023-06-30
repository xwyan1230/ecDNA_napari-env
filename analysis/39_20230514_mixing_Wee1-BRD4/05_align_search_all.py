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
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '6_Wee1-GFP-6hr_Ctrl-mCh-6hr'
batch = 2
total_fov = 40
dshape_factor = 0.0765
a_target = 1.1
n_col = 10

# topleft_target = [3313.627431155869, 1430.784821051334]  2_Ctrl-GFP-6hr_Ctrl-mCh-6hr [3333.627431155869, 1430.784821051334] 1.3
# topleft_target = [3753.8435086956088, 5124.393828360026]  5_Ctrl-GFP-6hr_Wee1-mCh-6hr [3733.8435086956088, 5094.393828360026] 1.1
# topleft_target = [5187.786656194921, 4282.431129459471]  6_Wee1-GFP-6hr_Ctrl-mCh-6hr [5177.786656194921, 4282.431129459471] 1.1
# topleft_target = [2615.5115503413313, 4631.470673742484]  9_Ctrl-GFP-6hr_BRD4-mCh-6hr [2635.5115503413313, 4641.470673742484] 0.5
# topleft_target = [3875.6283976329314, 2309.720277882853]  10_BRD4-GFP-6hr_Ctrl-mCh-6hr [3865.6283976329314, 2319.720277882853] 0.3

img_hoechst_stack = skio.imread("%s%s/%s_%s_%s_40pos_RAW_ch00.tif" % (data_dir, sample, exp, sample, batch), plugin="tifffile")

data = pd.DataFrame(columns=['sample', 'batch', 'fov', 'topleft_x', 'topleft_y', 'min_ratio'])

img_hoechst_Keyence = skio.imread("%s%s/Hoechst_cut.tif" % (data_dir, sample), plugin="tifffile")
img_hoechst_local = np.zeros_like(img_hoechst_Keyence)

fov = 0
print(fov)
img_hoechst = img_hoechst_stack[fov, :, :]
ori_shape = img_hoechst.shape
img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                interpolation=cv2.INTER_AREA)
img_hoechst_resize_rotate = imutils.rotate(img_hoechst_resize, angle=a_target)
s = img_hoechst_resize_rotate.shape[0]
border = 200
# topleft_target = [img_hoechst_Keyence.shape[0]-s-50-border, 50+border]
topleft_target = [img_hoechst_Keyence.shape[0]-s-50-border, 1535+1.05*s]
interval = 5
xrange = 200
yrange = 200
min_ratio = 1
topleft_ori1 = [topleft_target[0] - xrange/2, topleft_target[1] - yrange/2]
i = 0
for x in range(int(xrange/interval)):
    for y in range(int(yrange/interval)):
        i = i+1
        topleft = [topleft_ori1[0] + x * interval, topleft_ori1[1] + y * interval]
        img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
        img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft[1]):int(topleft[1] + s), int(topleft[0]):int(topleft[0] + s)]

        minus = img_hoechst_resize_rotate.astype(float) - img_hoechst_Keyence_cut.astype(float)
        minus[minus < 0] = 0
        min_ratio_temp = minus.sum() / img_hoechst_resize_rotate.sum()
        if min_ratio_temp < min_ratio:
            min_ratio = min_ratio_temp
            topleft_target = [topleft_ori1[0] + x * interval, topleft_ori1[1] + y * interval]
topleft0 = topleft_target
min_ratio0 = min_ratio
print(topleft0)
print(min_ratio0)

for fov in range(total_fov):
    print(fov)
    img_hoechst = img_hoechst_stack[fov, :, :]
    ori_shape = img_hoechst.shape
    img_hoechst_resize = cv2.resize(img_hoechst, dsize=(int(ori_shape[1] * dshape_factor), int(ori_shape[0] * dshape_factor)),
                                    interpolation=cv2.INTER_AREA)
    img_hoechst_resize_rotate = imutils.rotate(img_hoechst_resize, angle=a_target)
    s = img_hoechst_resize_rotate.shape[0]
    if fov == 0:
        topleft_target = topleft0
        data.loc[len(data.index)] = [sample, batch, 0, topleft0[0], topleft0[1], min_ratio0]
    else:
        if fov//n_col == 0:
            topleft_target = [int(topleft0[0]-((fov%n_col)*1.05*s)), int(topleft0[1]+((fov//n_col)*1.05*s))]
        else:
            ref_fov = fov-(2*(fov%n_col)+1)
            topleft_target = [data[data['fov'] == ref_fov]['topleft_x'].tolist()[0], int(data[data['fov'] == ref_fov]['topleft_y'].tolist()[0]+1.1*s)]
        print(topleft_target)
        min_ratio = 1
        topleft_ori1 = [topleft_target[0] - xrange/2, topleft_target[1] - yrange/2]
        i = 0
        for x in range(int(xrange / interval)):
            for y in range(int(yrange / interval)):
                i = i + 1
                topleft = [topleft_ori1[0] + x * interval, topleft_ori1[1] + y * interval]
                img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
                img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft[1]):int(topleft[1] + s),
                                     int(topleft[0]):int(topleft[0] + s)]

                minus = img_hoechst_resize_rotate.astype(float) - img_hoechst_Keyence_cut.astype(float)
                minus[minus < 0] = 0
                min_ratio_temp = minus.sum() / img_hoechst_resize_rotate.sum()
                if min_ratio_temp < min_ratio:
                    min_ratio = min_ratio_temp
                    topleft_target = [topleft_ori1[0] + x * interval, topleft_ori1[1] + y * interval]

        data.loc[len(data.index)] = [sample, batch, fov, topleft_target[0], topleft_target[1], min_ratio]
        print(min_ratio)
        print(topleft_target)

    img_hoechst_local = ima.image_paste_to(img_hoechst_local, img_hoechst_resize_rotate, [int(topleft_target[1]), int(topleft_target[0])])

    img_hoechst_Keyence_cut = img_hoechst_Keyence.copy()
    img_hoechst_Keyence_cut = img_hoechst_Keyence_cut[int(topleft_target[1]):int(topleft_target[1] + s), int(topleft_target[0]):int(topleft_target[0] + s)]
    data.to_csv('%s%s/%s_%s_alignment.txt' % (output_dir, sample, sample, batch), index=False, sep='\t')

    viewer = napari.Viewer()
    viewer.add_image(img_hoechst_Keyence_cut, blending='additive', colormap='blue', contrast_limits=[0, 45535])
    viewer.add_image(img_hoechst_resize_rotate, blending='additive', colormap='green', contrast_limits=[0, 45535])
    # viewer.add_image(minus, blending='additive')
    if not os.path.exists("%s%s/align_DNAFISH/" % (output_dir, sample)):
        os.makedirs("%s%s/align_DNAFISH/" % (output_dir, sample))
    plt.imsave("%s%s/align_DNAFISH/%s_fov%s_alignment_local.tiff" % (output_dir, sample, batch, fov), dis.blending(viewer))
    viewer.close()

viewer = napari.Viewer()
viewer.add_image(img_hoechst_Keyence, blending='additive', colormap='blue', contrast_limits=[0, 45535])
viewer.add_image(img_hoechst_local, blending='additive', colormap='green', contrast_limits=[0, 45535])
plt.imsave("%s%s/align_DNAFISH/%s_alignment_global.tiff" % (output_dir, sample, batch), dis.blending(viewer))
napari.run()