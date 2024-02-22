import napari
import shared.display as dis
import matplotlib.pyplot as plt
import seaborn as sns
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import seaborn as sns
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231026_analysis_ColoDMandHSR_Wee1/"

# hue_order = ['G1', 'S', 'G2', 'M']
# line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255), (0.85, 0.35, 0.25)]

sample = 'Colo320DM_Wee1-degrader'
conc = 'XY09'
fov = 1

df = pd.read_csv('%s/%s_%s_%s.txt' % (master_folder, sample, conc, fov), na_values=['.'], sep='\t')
df['ln_mean_int_green'] = np.log(df['green'])
df['ln_mean_int_red'] = np.log(df['red'])
df['ln_mean_int_mitosis'] = np.log(df['mitosis'])

df1 = df[df['ln_mean_int_mitosis'] >= 7].copy().reset_index(drop=True)

print(len(df))
print(len(df1))
print(len(df1[(df1['ln_mean_int_red'] < 7) & (df1['ln_mean_int_green'] < 6.5)]))
print(len(df1[(df1['ln_mean_int_red'] < 7) & (df1['ln_mean_int_green'] < 6.5)])/len(df1))
print(len(df1[(df1['ln_mean_int_red'] < 7) & (df1['ln_mean_int_green'] >= 6.5)]))
print(len(df1[(df1['ln_mean_int_red'] < 7) & (df1['ln_mean_int_green'] >= 6.5)])/len(df1))
print(len(df1[(df1['ln_mean_int_red'] >= 7) & (df1['ln_mean_int_green'] < 6.5)]))
print(len(df1[(df1['ln_mean_int_red'] >= 7) & (df1['ln_mean_int_green'] < 6.5)])/len(df1))
print(len(df1[(df1['ln_mean_int_red'] >= 7) & (df1['ln_mean_int_green'] >= 6.5)]))
print(len(df1[(df1['ln_mean_int_red'] >= 7) & (df1['ln_mean_int_green'] >= 6.5)])/len(df1))
print(len(df[df['ln_mean_int_mitosis'] > 9.5]))
print(len(df[df['ln_mean_int_mitosis'] > 9.5])/len(df1))

plt.subplots(figsize=(9, 6))
sns.distplot(df['ln_mean_int_mitosis'], bins=50, kde=False, hist_kws={'range':(5,11)})
plt.savefig('%s/pH3_%s_%s_%s.pdf' % (master_folder, sample, conc, fov))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df1, y='ln_mean_int_green', x='ln_mean_int_red', c=df1['ln_mean_int_mitosis'].tolist(), cmap='coolwarm', s=10, alpha=1, vmin=5, vmax=11)
plt.savefig('%s/%s_%s_%s.pdf' % (master_folder, sample, conc, fov))
plt.show()


