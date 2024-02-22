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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df_seq = pd.read_csv('%sseq_all.txt' % data_dir, na_values=['.'], sep='\t')

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]

temp = np.arange(0.0125, 1, 0.025)
column_lst = [str('%.3f' % i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)
data_heatmap1 = pd.DataFrame(columns=column_lst)
data_heatmap0 = pd.DataFrame(columns=column_lst)
data_heatmap_all = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_radial_summary_relativer.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')

    temp = []
    for j in column_lst:
        temp.append(np.log2(df[j][1] / (df[j][0] + 0.0001)))
    data_heatmap.loc[len(data_heatmap.index)] = temp
    data_heatmap1.loc[len(data_heatmap1.index)] = df.iloc[1]
    data_heatmap0.loc[len(data_heatmap0.index)] = df.iloc[0]
    data_heatmap_all.loc[len(data_heatmap_all.index)] = df.iloc[0]
    data_heatmap_all.loc[len(data_heatmap_all.index)] = df.iloc[1]

data_heatmap.index = df_seq['gene']
data_heatmap1.index = df_seq['gene']
data_heatmap0.index = df_seq['gene']
index_name = []
for i in range(len(df_seq)):
    index_name.append(df_seq['gene'][i])
    index_name.append('')
data_heatmap_all.index = index_name

"""plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_log2FC.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=2, vmin=-2, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_log2FC_fixscale.pdf' % output_dir)
plt.show()"""

"""plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap1, cbar=0, linewidths=2, vmax=data_heatmap1.values.max(), vmin=data_heatmap1.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer.pdf' % output_dir)
plt.show()"""

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap1, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_fixscale.pdf' % output_dir)
plt.show()

"""plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap0, cbar=0, linewidths=2, vmax=data_heatmap0.values.max(), vmin=data_heatmap0.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_ctrl.pdf' % output_dir)
plt.show()"""

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap0, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_ctrl_fixscale.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap_all, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_relativer_all_fixscale.pdf' % output_dir)
plt.show()

###############

temp = np.arange(98.75, 0, -2.5)
column_lst = [str('%.3f' % i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)
data_heatmap1 = pd.DataFrame(columns=column_lst)
data_heatmap0 = pd.DataFrame(columns=column_lst)
data_heatmap_all = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_radial_summary_absoluter.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')

    temp = []
    for j in column_lst:
        temp.append(np.log2(df[j][1] / (df[j][0] + 0.0001)))
    data_heatmap.loc[len(data_heatmap.index)] = temp
    data_heatmap1.loc[len(data_heatmap1.index)] = df.iloc[1]
    data_heatmap0.loc[len(data_heatmap0.index)] = df.iloc[0]
    data_heatmap_all.loc[len(data_heatmap_all.index)] = df.iloc[0]
    data_heatmap_all.loc[len(data_heatmap_all.index)] = df.iloc[1]

data_heatmap.index = df_seq['gene']
data_heatmap1.index = df_seq['gene']
data_heatmap0.index = df_seq['gene']
index_name = []
for i in range(len(df_seq)):
    index_name.append(df_seq['gene'][i])
    index_name.append('')
data_heatmap_all.index = index_name

"""plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_log2FC.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=2, vmin=-2, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_log2FC_fixscale.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap1, cbar=0, linewidths=2, vmax=data_heatmap1.values.max(), vmin=data_heatmap1.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter.pdf' % output_dir)
plt.show()"""

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap1, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_fixscale.pdf' % output_dir)
plt.show()

"""plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap0, cbar=0, linewidths=2, vmax=data_heatmap0.values.max(), vmin=data_heatmap0.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_ctrl.pdf' % output_dir)
plt.show()"""

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap0, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_ctrl_fixscale.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap_all, cbar=0, linewidths=2, vmax=2.0, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_seg_DNAFISH_absoluter_all_fixscale.pdf' % output_dir)
plt.show()



