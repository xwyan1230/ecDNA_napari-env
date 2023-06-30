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
df1 = pd.read_csv('%s%s_cellcycle_hoechstonly.txt' % (output_dir, sample), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['mean_int_green'])
df['ln_mean_int_red'] = np.log(df['mean_int_red'])
df['ln_mean_int_nuclear'] = np.log(df['mean_int_nuclear'])
df['ln_mean_int_MYC'] = np.log(df['mean_int_MYC'])
df1['ln_mean_int_green'] = np.log(df1['mean_int_green'])
df1['ln_mean_int_red'] = np.log(df1['mean_int_red'])
df1['ln_mean_int_nuclear'] = np.log(df1['mean_int_nuclear'])
df1['ln_mean_int_MYC'] = np.log(df1['mean_int_MYC'])

# df_M = df[df['intensity_stdev'] > 15000].copy()
# df_G1 = df[(df['ln_mean_int_red'] > 6) & (df['ln_mean_int_green'] < 5.5) & (df['intensity_stdev'] <= 15000)].copy()
# df_S = df[(df['ln_mean_int_red'] < 6) & (df['ln_mean_int_green'] > 5.5) & (df['intensity_stdev'] <= 15000)].copy()
# df_G2 = df[(df['ln_mean_int_red'] > 6) & (df['ln_mean_int_green'] > 5.5) & (df['intensity_stdev'] <= 15000)].copy()
# df_neg = df[(df['ln_mean_int_red'] < 6) & (df['ln_mean_int_green'] < 5.5) & (df['intensity_stdev'] <= 15000)].copy()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', color='red', s=5)
# plt.savefig('%s/%s_confirmedM.pdf' % (output_dir, sample))
plt.show()