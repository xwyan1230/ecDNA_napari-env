import nd2
import napari
import pandas as pd
import numpy as np
import shared.image as ima
import tifffile as tif
import skimage.io as skio

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240715_analysis_PCC_cellcycle/"
data_dir = '%sdata/' % master_folder

sample = '02_G1'
n_total = 60
for i in range(n_total):
    if i < 9:
        img_stack = nd2.imread('%s/%s/%s_00%s.nd2' % (data_dir, sample, sample, i+1))
    else:
        img_stack = nd2.imread('%s/%s/%s_0%s.nd2' % (data_dir, sample, sample, i+1))
    # print(img_stack.shape)

    # img1 = img_stack[:, 0, :, :]
    # img2 = img_stack[:, 1, :, :]
    img1 = img_stack[:, 0, :, :].max(axis=0)
    img2 = img_stack[:, 1, :, :].max(axis=0)

    viewer = napari.Viewer()
    viewer.add_image(img1, blending='additive', colormap='blue', contrast_limits=[0, 7000])
    viewer.add_image(img2, blending='additive', colormap='green', contrast_limits=[0, 10000])
    napari.run()
