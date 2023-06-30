import skimage.io as skio
import pandas as pd
from skimage.morphology import disk, dilation, medial_axis
from skimage.measure import label, regionprops
import shared.image as ima
import numpy as np
import shared.dataframe as dat
import shared.math as mat
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_DNAFISH_ctrl/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

exp = '20230409_DM_plate_DNAFISH_control'
df = pd.read_csv('%s/summary.txt' % data_dir1, na_values=['.'], sep='\t')

x = 'mean_int_DNAFISH'
y = 'averageD_0.5'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 1)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH'
y = 'n_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 10)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH'
y = 'averageD_mean'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 1)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH'
y = 'averageD_std'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 0.4)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH'
y = 'total_area_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 2000)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH'
y = 'total_area_ratio_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
plt.xlim([0, 10000])
plt.ylim([0, max(df[y].max(), 0.15)])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()