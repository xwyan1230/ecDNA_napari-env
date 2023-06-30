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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sdata/" % master_folder

sample = 'DM_324pos_merge'
df = pd.read_csv('%s%s_cellcycle_updated.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])
df['ln_mean_int_nuclear'] = np.log(df['mean_int_nuclear'])
df['ln_mean_int_MYC'] = np.log(df['mean_int_MYC'])

df['ln_mean_int_red_cal'] = df['ln_mean_int_red'] - 0.14 * df['ln_mean_int_green']
df['ln_mean_int_green_cal'] = df['ln_mean_int_green'] - 0.2 * df['ln_mean_int_red']
df_G1 = df[(df['ln_mean_int_red_cal'] > 5) & (df['ln_mean_int_green_cal'] < 5) & (df['cellcycle_updated'] != 'M')].copy()
df_S = df[(df['ln_mean_int_red_cal'] < 5) & (df['ln_mean_int_green_cal'] > 5) & (df['cellcycle_updated'] != 'M')].copy()
df_G2 = df[(df['ln_mean_int_red_cal'] > 5) & (df['ln_mean_int_green_cal'] > 5) & (df['cellcycle_updated'] != 'M')].copy()
df_neg = df[(df['ln_mean_int_red_cal'] < 5) & (df['ln_mean_int_green_cal'] < 5) & (df['cellcycle_updated'] != 'M')].copy()

cellcycle = []
for i in range(len(df)):
    if df['cellcycle_updated'][i] == 'M':
        cellcycle.append('M')
    elif (df['ln_mean_int_red_cal'][i] > 5) & (df['ln_mean_int_green_cal'][i] < 5):
        cellcycle.append('G1')
    elif (df['ln_mean_int_red_cal'][i] < 5) & (df['ln_mean_int_green_cal'][i] > 5):
        cellcycle.append('S')
    elif (df['ln_mean_int_red_cal'][i] > 5) & (df['ln_mean_int_green_cal'][i] > 5):
        cellcycle.append('G2')
    elif (df['ln_mean_int_red_cal'][i] < 5) & (df['ln_mean_int_green_cal'][i] < 5):
        cellcycle.append('neg')
    else:
        cellcycle.append('NA')
df['cellcycle'] = cellcycle
df.to_csv('%s%s_cellcycle_updated1.txt' % (output_dir, sample), index=False, sep='\t')
df_temp = df[df['cellcycle'].isin(['G1', 'S', 'G2', 'M'])].copy().reset_index(drop=True)
df = df_temp

print(len(df_G1))
print(len(df_S))
print(len(df_G2))
print(len(df[df['cellcycle_updated'] == 'M']))

plt.subplots(figsize=(6, 6))
sns.scatterplot(data=df, y='ln_mean_int_green_cal', x='ln_mean_int_red_cal', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
sns.scatterplot(data=df[df['cellcycle_updated'] == 'M'], y='ln_mean_int_green_cal', x='ln_mean_int_red_cal', color='red', s=5)
# plt.savefig('%s/%s_confirmedM.pdf' % (output_dir, sample))
plt.show()