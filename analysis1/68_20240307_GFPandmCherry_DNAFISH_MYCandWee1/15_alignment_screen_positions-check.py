import skimage.io as skio
import napari
import cv2
import shared.display as dis
import matplotlib.pyplot as plt
import numpy as np
import imutils
import pandas as pd
import shared.dataframe as dat
import shared.image as ima
import tifffile as tif
import os

# INPUT PARAMETERS
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240307_analysis_GFPandmCherry_DNAFISH/"
sample = 'F4'
afterFISH_hoechst_pixel = 766
DNAFISH_hoechst_fov_pixel = 180

# NO NEED TO CHANGE
data_dir = "%sdata/" % master_folder
data_dir1 = "%sprocessed/" % master_folder
align_dir = "%salign/" % master_folder
output_dir = "%sprocessed/" % master_folder

# DO NOT CHANGE
dshape_factor = DNAFISH_hoechst_fov_pixel * 1.0 / afterFISH_hoechst_pixel
border = 100

align_topleft = pd.read_csv('%s/align_screen_topleft_locations_%s.txt' % (align_dir, sample), na_values=['.'], sep='\t')
align_topleft['topleft_target'] = [dat.str_to_float(x) for x in align_topleft['topleft_target']]
img_after_hoechst = skio.imread("%s/%s/%s_afterFISH_hoechst_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")
img_screen_cut = skio.imread("%s/%s/%s_screen_cut.tif" % (data_dir1, sample, sample), plugin="tifffile")

viewer = napari.Viewer()
viewer.add_image(img_after_hoechst, blending='additive', colormap='blue', contrast_limits=[0, 65535])
viewer.add_image(img_screen_cut, blending='additive', contrast_limits=[0, img_screen_cut.max()])
points = [[int(align_topleft['topleft_target'][x][0]), int(align_topleft['topleft_target'][x][1])] for x in range(len(align_topleft))]
points_layer = viewer.add_points(points, size=10)
napari.run()
