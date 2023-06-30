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

hue_order = ['prophase', 'metaphase', 'anaphase', 'telophase']

df = pd.read_csv('%s%s_cellcycle.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df_M = pd.read_csv('%s%s_mitosis.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df_M1 = df_M[df_M['mitosis'].isin(hue_order)].copy().reset_index(drop=True)

mitosis = []
cellcycle_updated = []
for i in range(len(df)):
    if df['label'][i] in df_M['label'].tolist():
        if df['label'][i] in df_M1['label'].tolist():
            mitosis.append(df_M1[df_M1['label'] == df['label'][i]]['mitosis'].tolist()[0])
            cellcycle_updated.append('M')
        else:
            mitosis.append('NA')
            cellcycle_updated.append('debris')
    else:
        mitosis.append(df['cellcycle'][i])
        cellcycle_updated.append(df['cellcycle'][i])
df['mitosis'] = mitosis
df['cellcycle_updated'] = cellcycle_updated
df.to_csv('%s%s_cellcycle_updated.txt' % (output_dir, sample), index=False, sep='\t')
print("DONE!")

