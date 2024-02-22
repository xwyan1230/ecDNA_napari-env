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

df_seq = pd.read_csv('%sseq.txt' % data_dir, na_values=['.'], sep='\t')

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]

temp = np.arange(0.0125, 1, 0.025)
column_lst = [str('%.3f' % i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH', 'radial_curve_edge_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

    data_heatmap_temp = pd.DataFrame(columns=column_lst)
    for s in hue_order:
        data_sample = df[df['group'] == s].copy().reset_index(drop=True)
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]
        data_heatmap_temp.loc[len(data_heatmap_temp.index)] = data_radial.mean()

    temp = []
    for j in column_lst:
        temp.append(np.log2(data_heatmap_temp[j][1] / (data_heatmap_temp[j][0] + 0.0001)))
    data_heatmap.loc[len(data_heatmap.index)] = temp

data_heatmap.index = df_seq['gene']

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_DNAFISH.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=0.2, vmin=-0.2, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_DNAFISH_fixscale.pdf' % output_dir)
plt.show()

###############

temp = np.arange(0.0125, 1, 0.025)
column_lst = [str('%.3f' % i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH', 'radial_curve_edge_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

    for s in hue_order:
        data_sample = df[df['group'] == s].copy().reset_index(drop=True)
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]
        data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

index_name = []
for i in range(len(df_seq)):
    index_name.append(df_seq['gene'][i])
    index_name.append('')
data_heatmap.index = index_name

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_all_DNAFISH.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_all_DNAFISH_fixscale.pdf' % output_dir)
plt.show()

###############

temp = np.arange(98.75, 0, -2.5)
column_lst = [str(i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH', 'radial_curve_edge_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

    data_heatmap_temp = pd.DataFrame(columns=column_lst)
    for s in hue_order:
        data_sample = df[df['group'] == s].copy().reset_index(drop=True)
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [data_sample['radial_curve_edge_DNAFISH'][i][-(j+1)] for i in range(len(data_sample))]
        data_heatmap_temp.loc[len(data_heatmap_temp.index)] = data_radial.mean()

    temp = []
    for j in column_lst:
        temp.append(np.log2(data_heatmap_temp[j][1] / (data_heatmap_temp[j][0] + 0.0001)))
    data_heatmap.loc[len(data_heatmap.index)] = temp

data_heatmap.index = df_seq['gene']

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_edge_DNAFISH.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=0.2, vmin=-0.2, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_edge_DNAFISH_fixscale.pdf' % output_dir)
plt.show()

#################

temp = np.arange(98.75, 0, -2.5)
column_lst = [str(i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)

for i in range(len(df_seq)):
    sample = df_seq['location'][i]
    df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH', 'radial_curve_edge_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

    for s in hue_order:
        data_sample = df[df['group'] == s].copy().reset_index(drop=True)
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [data_sample['radial_curve_edge_DNAFISH'][i][-(j+1)] for i in range(len(data_sample))]
        data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

index_name = []
for i in range(len(df_seq)):
    index_name.append(df_seq['gene'][i])
    index_name.append('')
data_heatmap.index = index_name

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_edge_all_DNAFISH.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(16, 16))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_line/heatmap_radial_edge_all_DNAFISH_fixscale.pdf' % output_dir)
plt.show()

