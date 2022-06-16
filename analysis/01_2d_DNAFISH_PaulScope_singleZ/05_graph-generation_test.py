import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import skimage.io as skio
import napari

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
raw_folder = '01_raw'
seg_folder = '02_seg'
total_fov = 16
total_z = 18

for fov in range(total_fov):
    im_z_stack_DNAFISH = skio.imread("%s%s/%s/%s_RAW_ch01_fov%s.tif" % (master_folder, sample, raw_folder, sample, fov),
                                     plugin="tifffile")
    max_DNAFISH = im_z_stack_DNAFISH.max(axis=0)
    std_DNAFISH = im_z_stack_DNAFISH.std(axis=0)

    viewer = napari.view_image(std_DNAFISH)
    napari.run()
