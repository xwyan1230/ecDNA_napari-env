import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230514_analysis_mixing_Wee1-BRD4/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230512_mixing_Wee1-BRD4_6hr'
sample = '10_BRD4-GFP-6hr_Ctrl-mCh-6hr'

df = pd.read_csv('%s%s/%s_n4.txt' % (data_dir1, sample, sample), na_values=['.'], sep='\t')
df_label = pd.read_csv('%s%s/alignment_label.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
df_keyence = pd.read_csv('%s%s/analysis_keyence.txt' % (data_dir1, sample), na_values=['.'], sep='\t')

df['keyence_index'] = df_label['index']
df_sort = df[df['keyence_index'] != 0].copy().reset_index(drop=True)

int_GFP = []
int_mCherry = []
group = []

for i in range(len(df_sort)):
    print('%s: %s' % (i, df_sort['keyence_index'][i]))
    df_temp = df_keyence[df_keyence['label'] == df_sort['keyence_index'][i]]
    int_GFP.append(df_temp['mean_int_GFP'].tolist()[0])
    int_mCherry.append(df_temp['mean_int_mCherry'].tolist()[0])
    group.append(df_temp['group'].tolist()[0])

df_sort['mean_int_GFP'] = int_GFP
df_sort['mean_int_mCherry'] = int_mCherry
df_sort['group'] = group

df_sort.to_csv('%s%s/%s_summary.txt' % (output_dir, sample, sample), index=False, sep='\t')
print("DONE!")