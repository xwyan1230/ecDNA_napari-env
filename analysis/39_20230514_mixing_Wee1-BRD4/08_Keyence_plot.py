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
df = pd.read_csv('%s%s/analysis_keyence.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
df['ln_mean_int_GFP'] = np.log(df['mean_int_GFP'])
df['ln_mean_int_mCherry'] = np.log(df['mean_int_mCherry'])

df_mCherry = df[(df['ln_mean_int_mCherry'] > 7.5) & (df['ln_mean_int_GFP'] < 8.4)].copy().reset_index(drop=True)
df_GFP = df[(df['ln_mean_int_mCherry'] < 7) & (df['ln_mean_int_GFP'] > 8.4)].copy().reset_index(drop=True)

group = []
for i in range(len(df)):
    if (df['ln_mean_int_mCherry'][i] > 7.5) & (df['ln_mean_int_GFP'][i] < 8.4):
        group.append('mCherry')
    elif (df['ln_mean_int_mCherry'][i] < 7) & (df['ln_mean_int_GFP'][i] > 8.4):
        group.append('GFP')
    else:
        group.append('NA')

df['group'] = group
df.to_csv('%s%s/analysis_keyence.txt' % (output_dir, sample), index=False, sep='\t')

print(len(df))
print(len(df_mCherry))
print(len(df_GFP))

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='ln_mean_int_GFP', x='ln_mean_int_mCherry', color=(0.30, 0.30, 0.30), s=5, alpha=0.3)
sns.scatterplot(data=df_mCherry, y='ln_mean_int_GFP', x='ln_mean_int_mCherry', color='red', s=5, alpha=0.5)
sns.scatterplot(data=df_GFP, y='ln_mean_int_GFP', x='ln_mean_int_mCherry', color='green', s=5, alpha=0.5)
plt.savefig('%s%s/group_keyence_total%s_GFP%s_mCh%s.pdf' % (output_dir, sample, len(df), len(df_GFP), len(df_mCherry)))
plt.show()