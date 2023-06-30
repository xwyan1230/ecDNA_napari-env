import napari
import shared.display as dis
import tifffile as tif
import matplotlib.pyplot as plt
import numpy as np
import cv2
import math
import shared.segmentation as seg
from shared.sinaplot import sinaplot
import shared.objects as obj
import pandas as pd
import seaborn as sns
import shared.image as ima
from skimage.measure import label, regionprops_table, regionprops
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk
import os
import skimage.io as skio

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'

hue_order = ['G1', 'S', 'G2', 'prophase', 'metaphase', 'anaphase', 'telophase']

df = pd.read_csv('%s%s_cellcycle_updated1.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df['total_int_MYC'] = df['area_nuclear']*df['mean_int_MYC']

df_cc = df[df['mitosis'].isin(hue_order)]

# sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
feature = 'total_int_MYC'
sinaplot(data=df, x='mitosis', y=feature, order=hue_order, violin=False, scale='area', point_size=2)
plt.ylim([0, 3E8])
plt.savefig('%s/%s_%s_cellcycle1.pdf' % (output_dir, sample, feature))
plt.show()

plt.subplots(figsize=(9, 6))
feature = 'total_int_MYC'
sinaplot(data=df, x='cellcycle_updated', y=feature, order=['G1', 'S', 'G2', 'M'], violin=False, scale='area', point_size=2)
plt.ylim([0, 3E8])
plt.savefig('%s/%s_%s2.pdf' % (output_dir, sample, feature))
plt.show()