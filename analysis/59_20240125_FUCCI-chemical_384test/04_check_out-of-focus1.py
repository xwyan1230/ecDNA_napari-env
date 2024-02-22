import skimage.io as skio
from skimage.filters import threshold_otsu
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240125_analysis_FUCCI-chemical_384test/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

# set parameters
otsu_factor = 1.1
circ_threshold = 0.5
min_size = 300
max_size = 1500

folder = '20x_3x3_2500'
sample = 'XY02'
n_img = 9

data = pd.DataFrame()
fov_lst = []
hoechst_lst = []
h23_lst = []
azaleaB5_lst = []
emiRFP670_lst = []


def img_crop(img, i):
    if i == 0:
        out = img[0:1440, 0:1920]
    elif (i == 1) | (i == 2):
        out = img[0:1440, 576:1920]
    elif (i == 3) | (i == 6):
        out = img[432:1440, 0:1920]
    elif (i == 4) | (i == 5):
        out = img[432:1440, 0:1344]
    elif (i == 7) | (i == 8):
        out = img[432:1440, 576:1920]
    return out


for i in range(n_img):
    print(i)
    file_name = 'Image_%s_0000%s' % (sample, i+1)
    img_hoechst = img_crop(skio.imread("%s%s/%s/%s_CH1.tif" % (data_dir, folder, sample, file_name), plugin="tifffile")[:, :, 2], i)
    thresh_val = threshold_otsu(img_hoechst)
    print(thresh_val)

