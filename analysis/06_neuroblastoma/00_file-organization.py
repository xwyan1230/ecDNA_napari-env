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
import pandas as pd
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220826_neuroblastoma/"
group = 'B'
group_folder = '%s%s/' % (master_folder, group)
save_path = master_folder

multi_imgs = [x for x in os.listdir(group_folder)]
if '.DS_Store' in multi_imgs:
    multi_imgs.remove('.DS_Store')

for i in multi_imgs:
    sample = i.replace('-', '_').split('.')[0].split('_')[0]
    color_string = i
    color = 'color'
    if 'r' in color_string:
        color = 'r'
    elif 'g' in color_string:
        color = 'g'
    elif 'b' in color_string:
        color = 'b'
    if color != 'color':
        os.rename('%s%s' % (group_folder, i), '%s%s_%s.BMP' % (group_folder, sample, color))

print("DONE!")
