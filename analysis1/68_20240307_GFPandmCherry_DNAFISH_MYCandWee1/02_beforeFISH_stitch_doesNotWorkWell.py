import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import m2stitch
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 2000
# threshold = 3000
n_img = 88

name = pd.read_csv('%s/name.txt' % data_dir, na_values=['.'], sep='\t')

samples_lst = name['sample'].tolist()
wells_lst = name['well'].tolist()
beforeFISH_lst = name['beforeFISH'].tolist()
afterFISH_lst = name['afterFISH'].tolist()

# https://m2stitch.readthedocs.io/en/latest/usage.html

rows = [1] * 8 + [2] * 8 + [3] * 8 + [4] * 8 + [5] * 8 + [6] * 8 + [7] * 8 + [8] * 8 + [9] * 8 + [10] * 8 + [11] * 8
cols = list(np.arange(1, 9, 1)) + list(np.arange(8, 0, -1)) + list(np.arange(1, 9, 1)) + list(np.arange(8, 0, -1)) + \
       list(np.arange(1, 9, 1)) + list(np.arange(8, 0, -1)) + list(np.arange(1, 9, 1)) + list(np.arange(8, 0, -1)) + \
       list(np.arange(1, 9, 1)) + list(np.arange(8, 0, -1)) + list(np.arange(1, 9, 1))

for k in range(len(samples_lst)):
    sample = samples_lst[k]
    well = wells_lst[k]
    beforeFISH = beforeFISH_lst[k]
    print(sample)
    print(well)

    img_3d_green = skio.imread("%s/beforeFISH/%s/Image_%s_00001_CH2.tif" % (data_dir, beforeFISH, beforeFISH), plugin="tifffile")[
        np.newaxis, :, :, 1]
    img_3d_red = skio.imread("%s/beforeFISH/%s/Image_%s_00001_CH3.tif" % (data_dir, beforeFISH, beforeFISH), plugin="tifffile")[
        np.newaxis, :, :, 0]

    for i in range(n_img-1):
        print(i)
        if i < 8:
            file_name = 'Image_%s_0000%s' % (beforeFISH, i+2)
        else:
            file_name = 'Image_%s_000%s' % (beforeFISH, i+2)
        img_green = skio.imread("%s/beforeFISH/%s/%s_CH2.tif" % (data_dir, beforeFISH, file_name), plugin="tifffile")[:, :, 1]
        img_red = skio.imread("%s/beforeFISH/%s/%s_CH3.tif" % (data_dir, beforeFISH, file_name), plugin="tifffile")[:, :, 0]
        img_3d_green = np.vstack([img_3d_green, img_green[np.newaxis, ...]])
        img_3d_red = np.vstack([img_3d_red, img_red[np.newaxis, ...]])
        print(img_3d_green.shape)

    result_df, _ = m2stitch.stitch_images(img_3d_green, rows, cols, row_col_transpose=False)

    result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
    result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

    size_y = img_3d_green.shape[1]
    size_x = img_3d_green.shape[2]

    stitched_image_size = (
        result_df["y_pos2"].max() + size_y,
        result_df["x_pos2"].max() + size_x,
    )
    stitched_image = np.zeros_like(img_3d_green, shape=stitched_image_size)
    for i, row in result_df.iterrows():
        stitched_image[
        row["y_pos2"]: row["y_pos2"] + size_y,
        row["x_pos2"]: row["x_pos2"] + size_x,
        ] = img_3d_green[i]

    viewer = napari.Viewer()
    viewer.add_image(stitched_image, blending='additive', colormap='green', contrast_limits=[0, 65535])
    napari.run()

print("DONE!")

