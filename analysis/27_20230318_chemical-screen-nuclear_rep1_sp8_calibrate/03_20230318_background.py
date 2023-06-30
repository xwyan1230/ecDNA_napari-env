import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, dilation
import shared.dataframe as dat
from skimage import segmentation
import shared.math as mat
import math
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

plate = 'DM_6hr'
cal = pd.read_csv('%s%s_calibration.txt' % (data_dir, plate), na_values=['.'], sep='\t')

samples = cal['sample'].tolist()

bg_hoechst_llst = []
bg_MYC_llst = []
bg_hoechst_mean_lst = []
bg_MYC_mean_lst = []

for s in range(len(samples)):
    sample = samples[s]

    # file_name = '20230318_HSR_6hr_rep1_calibration_%s_RAW' % (sample)
    file_name = '20230319_DM_6hr_rep1_calibrate_%s_RAW' % (sample)
    img_hoechst_stack = skio.imread("%s%s/%s_ch01.tif" % (data_dir, plate, file_name), plugin="tifffile")
    img_MYC_stack = skio.imread("%s%s/%s_ch00.tif" % (data_dir, plate, file_name), plugin="tifffile")

    for fov in range(4):
        bg_hoechst_lst = []
        bg_MYC_lst = []

        img_hoechst = img_hoechst_stack[:, :, fov]
        img_MYC = img_MYC_stack[:, :, fov]
        # img_MYC = np.concatenate([np.zeros(shape=[5, 2048]), img_MYC], axis=0)[:2048, :2048]

        viewer = napari.Viewer()
        viewer.add_image(img_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
        viewer.add_image(img_MYC, blending='additive', colormap='magenta', contrast_limits=[0, 65535])
        shapes = viewer.add_shapes(name='Shapes', ndim=2)
        napari.run()

        poly_data = shapes.data[0]
        shapes_layer = viewer.layers['Shapes']
        top, left = np.floor(np.min(poly_data, axis=0))
        bottom, right = np.ceil(np.max(poly_data, axis=0))
        top, bottom = np.clip((top, bottom), 0, img_hoechst.shape[0] - 1).astype(int)
        left, right = np.clip((left, right), 0, img_hoechst.shape[1] - 1).astype(int)
        output_shape = (bottom - top + 1, right - left + 1)
        # generate sub_masks and sub_channels
        sub_masks = shapes_layer._data_view.to_masks(mask_shape=output_shape, offset=(top, left))[0]
        img_hoechst_mask = img_hoechst[top:bottom + 1, left:right + 1] * sub_masks
        img_MYC_mask = img_MYC[top:bottom + 1, left:right + 1] * sub_masks
        bg_hoechst = regionprops(label(sub_masks), img_hoechst_mask)[0].intensity_mean
        bg_MYC = regionprops(label(sub_masks), img_MYC_mask)[0].intensity_mean

        bg_hoechst_lst.append(bg_hoechst)
        bg_MYC_lst.append(bg_MYC)

    bg_hoechst_llst.append(bg_hoechst_lst)
    bg_hoechst_mean_lst.append(np.mean(bg_hoechst_lst))
    bg_MYC_llst.append(bg_MYC_lst)
    bg_MYC_mean_lst.append(np.mean(bg_MYC_lst))

cal['bg_hoechst'] = bg_hoechst_llst
cal['bg_MYC'] = bg_MYC_llst
cal['bg_hoechst_mean'] = bg_hoechst_mean_lst
cal['bg_MYC_mean'] = bg_MYC_mean_lst

cal.to_csv('%s%s_calibration.txt' % (output_dir, plate), index=False, sep='\t')

print("DONE!")