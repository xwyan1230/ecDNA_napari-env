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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_6hr/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df_ctrl = pd.read_csv('%s/summary_ctrl.txt' % data_dir1, na_values=['.'], sep='\t')
df = pd.read_csv('%s/summary.txt' % data_dir1, na_values=['.'], sep='\t')

offset = 250
df_ctrl['mean_int_DNAFISH_cal'] = df_ctrl['mean_int_DNAFISH']+offset
df['mean_int_DNAFISH_cal'] = df['mean_int_DNAFISH']

df_temp = df_ctrl[df_ctrl['n'] > 50].copy().reset_index(drop=True)
df_ctrl = df_temp
df_temp = df[df['n'] > 50].copy().reset_index(drop=True)
df = df_temp

ctrls = ['C3', 'C10', 'D6', 'F3', 'F10']
df_WT = pd.DataFrame()
for i in range(len(ctrls)):
    df_temp = df[df['sample'] == ctrls[i]].copy().reset_index(drop=True)
    df_WT = pd.concat([df_WT, df_temp], axis=0).reset_index(drop=True)

x = 'mean_int_DNAFISH_cal'
y = 'averageD_0.5'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 1, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH_cal'
y = 'n_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 10, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH_cal'
y = 'averageD_mean'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 1, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH_cal'
y = 'averageD_std'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 0.4, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH_cal'
y = 'total_area_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 2000, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()

x = 'mean_int_DNAFISH_cal'
y = 'total_area_ratio_ecDNA'
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_ctrl, x=x, y=y, color=(0.8, 0.8, 0.8), s=20)
sns.scatterplot(data=df, x=x, y=y, color=(0.30, 0.30, 0.30), s=20)
for i in range(len(df)):
    plt.text(x=df[x][i]+0.001, y=df[y][i]+0.001, s=df['sample'][i], size=5)
sns.scatterplot(data=df_WT, x=x, y=y, color=(0.85, 0.35, 0.25), s=20)
plt.xlim([0, max(df[x].max(), 10000, df_ctrl[x].max())])
plt.ylim([0, max(df[y].max(), 0.15, df_ctrl[y].max())])
plt.savefig('%s/%s_vs_%s.pdf' % (output_dir, x, y))
plt.show()