import skimage.io as skio
import napari
import imutils
import shared.image as ima
import shared.objects as obj
from skimage.morphology import disk, dilation, medial_axis
import shared.dataframe as dat
import math
import shared.display as dis
from skimage.morphology import disk, dilation
import tifffile as tif
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops
import numpy as np
import pandas as pd
# from napari_animation import Animation
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sdata/raw/" % master_folder
data_dir1 = "%sdata/seg_tif/" % master_folder
data_dir2 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'
# prefix = '20230614_acidFISH_lamin_ColoDM_DM'
prefix = '20230614_acidFISH_lamin_ColoHSR_HSR'

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25), (154/255, 205/255, 50/255)]
name_order = ['laminB1 IF', 'MYC DNA FISH', 'hoechst stain']
hue_order = ['radial_curve_laminB', 'radial_curve_DNAFISH', 'radial_curve_nuclear']
hue_order1 = ['radial_curve_edge_laminB', 'radial_curve_edge_DNAFISH', 'radial_curve_edge_nuclear']
"""name_order = ['laminB1 IF', 'MYC DNA FISH']
hue_order = ['radial_curve_laminB', 'radial_curve_DNAFISH']
hue_order1 = ['radial_curve_edge_laminB', 'radial_curve_edge_DNAFISH']"""

df = pd.read_csv('%s%s/%s_n4.txt' % (data_dir2, sample, sample), na_values=['.'], sep='\t')

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_edge_nuclear', 'radial_curve_edge_DNAFISH',
           'percentage_area_curve_ecDNA', 'n_ecDNA_lst', 'radial_curve_laminB', 'radial_curve_edge_laminB',
           'total_area_DNAFISH_lst', 'total_area_ratio_DNAFISH_lst']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# radial heatmap
print("Plotting total radial curve...")
temp = np.arange(0.0125, 1, 0.025)
column_lst = [str('%.3f' % i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)
for s in range(len(hue_order)):
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [df[hue_order[s]][i][j] for i in range(len(df))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = name_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial curve
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(df[hue_order[k]].tolist())
    for i in range(len(df)):
        plt.plot(x, df[hue_order[k]][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (name_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0.2, 2.0])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
number_nuclear = len(df)
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(df[hue_order[1]].tolist())
for i in range(len(df)):
    plt.plot(x, df[hue_order[1]][i], alpha=0.01, color=line_colors[1])
plt.plot(x, mean_curve3, color=line_colors[1], label='%s, n=%s' % (name_order[1], number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISHonly.pdf' % (output_dir, sample, sample))
plt.close()

# radial heatmap
print("Plotting total radial curve...")
temp = np.arange(59.25, 0, -1.5)
column_lst = [str(i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)
for s in range(len(hue_order)):
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [df[hue_order1[s]][i][-(j+1)] for i in range(len(df))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = name_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_edge_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_edge_DNAFISH_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial curve
x = np.arange(59.25, 0, -1.5)
x_label = 'r to edge (reverse)'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(df[hue_order1[k]].tolist())
    for i in range(len(df)):
        plt.plot(x, df[hue_order1[k]][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (name_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0.2, 2.0])
plt.ylabel('radial_curve_edge')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_edge_DNAFISH.pdf' % (output_dir, sample, sample))
plt.close()

x = np.arange(59.25, 0, -1.5)
x_label = 'r to edge (reverse)'
plt.subplots(figsize=(12, 9))
number_nuclear = len(df)
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(df[hue_order1[1]].tolist())
for i in range(len(df)):
    plt.plot(x, df[hue_order1[1]][i], alpha=0.01, color=line_colors[1])
plt.plot(x, mean_curve3, color=line_colors[1], label='%s, n=%s' % (name_order[1], number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve_edge')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_edge_DNAFISHonly.pdf' % (output_dir, sample, sample))
plt.close()