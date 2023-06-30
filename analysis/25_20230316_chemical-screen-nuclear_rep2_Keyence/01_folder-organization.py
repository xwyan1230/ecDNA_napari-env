import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import shared.objects as obj
import shared.image as ima
from skimage import segmentation
import numpy as np
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import shutil
import pandas as pd
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230316_analysis_chemical-screen-nuclear_rep2_Keyence/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

plate = 'HSR_2hr'

if not os.path.exists("%s%s_mod/raw_img/" % (output_dir, plate)):
    os.makedirs("%s%s_mod/raw_img/" % (output_dir, plate))

if not os.path.exists("%s%s_mod/color_img/" % (output_dir, plate)):
    os.makedirs("%s%s_mod/color_img/" % (output_dir, plate))

multi_imgs_dir = [x[0] for x in os.walk('%s%s' % (data_dir, plate))]
multi_imgs_dir.remove('%s%s' % (data_dir, plate))
for i in multi_imgs_dir:
    img = i.split('/')[-1]
    shutil.move('%s/Image_%s_CH1.tif' % (i, img), '%s%s_mod/raw_img/Image_%s_CH1.tif' % (output_dir, plate, img))
    shutil.move('%s/Image_%s_CH3.tif' % (i, img), '%s%s_mod/raw_img/Image_%s_CH3.tif' % (output_dir, plate, img))
    shutil.move('%s/Image_%s_Overlay.tif' % (i, img), '%s%s_mod/color_img/Image_%s_Overlay.tif' % (output_dir, plate, img))

print("DONE!")
