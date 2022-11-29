import shared.segmentation as seg
from skimage.measure import regionprops, label
from skimage.filters import threshold_otsu, threshold_local, sobel
from skimage.morphology import extrema, binary_dilation, binary_erosion, disk, medial_axis
import shared.objects as obj
import math
import shared.image as ima
import shared.dataframe as dat
from skimage import segmentation
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
import seaborn as sns
import skimage.io as skio
import shared.math as mat
import shared.display as dis
import pandas as pd
import tifffile as tif
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221017_periphery-localization_analysis/20220826_neuroblastoma/E/"
sample_lst = ['AO',
              'BD', 'BQ', 'CK', 'CO', 'CR', 'DA', 'DT', 'FD', 'FM', 'GP',
              'GR', 'HD', 'HZ', 'IK', 'IR', 'JF', 'JZ', 'KK', 'KW', 'LM',
              'LT', 'LX', 'LY', 'ME', 'MK', 'MM', 'NB', 'NP', 'NQ']
sample_lst = ['DT', 'JF', 'NQ', 'GR', 'LX', 'DA', 'MK', 'LT', 'NB', 'JZ', 'LY', 'MM', 'GP', 'LM', 'CO', 'BQ', 'NP',
              'HD', 'BD', 'CR', 'IK', 'AO', 'CK', 'ME', 'FD', 'KK', 'KW', 'FM', 'HZ']
sample_lst = ['HZ', 'FM', 'KW', 'FD', 'ME', 'CK', 'AO', 'IK', 'CR', 'BD', 'HD', 'NP', 'BQ', 'CO', 'LM', 'GP',
              'MM', 'LY', 'JZ', 'NB', 'LT', 'MK', 'DA', 'LX', 'GR', 'NQ', 'JF', 'DT']
sample_lst = ['HZ', 'AO', 'CR', 'FD', 'KW', 'CK', 'IK', 'BQ', 'CO', 'LM', 'HD', 'BD', 'NP', 'GP',
              'MM', 'JZ', 'NB', 'LY', 'LT', 'MK', 'DA']
sample_lst = ['CR', 'AO', 'KW', 'BQ', 'IK', 'NP', 'HD', 'BD', 'CO', 'LM', 'CK', 'GP', 'LX', 'DA', 'MK']
sample_lst = ['DA', 'MK', 'NP', 'IK', 'HD', 'LX', 'KW', 'AO', 'CR', 'BQ', 'CO', 'GP', 'CK', 'LM', 'BD' ]
# n [13, 35, 19, 28, 11, 16, 37, 31, 21, 14, 15, 11, 15, 52, 17, 28, 25, 8, 15, 11, 29, 40, 8, 35, 11, 11, 10, 20, 25]

group = 'E'
save_path = master_folder
save_folder = master_folder
version = 1

feature = ['mean_int_ind_ecDNA', 'total_int_ind_ecDNA', 'area_ind_ecDNA', 'area_ratio_ind_ecDNA', 'g',
           'radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized', 'radial_curve_mask',
           'angle_curve_nuclear', 'angle_curve_DNAFISH',
           'percentage_area_curve_ecDNA', 'percentage_area_ratio_curve_ecDNA', 'percentage_int_curve_ecDNA',
           'cum_area_ind_ecDNA', 'cum_area_ratio_ind_ecDNA', 'cum_int_ind_ecDNA', 'cum_area_ind_ecDNA_filled',
           'cum_area_ratio_ind_ecDNA_filled', 'cum_int_ind_ecDNA_filled']

# LOAD FILE
data = pd.read_csv('%s%s_v%s.txt' % (master_folder, group, version), na_values=['.'], sep='\t')

for f in feature:
    data[f] = [dat.str_to_float(data[f][i]) for i in range(len(data))]

# heatmap
data_heatmap = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
data_heatmap_mean = pd.DataFrame(columns=['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
sample_updated = []
sample_n = []
for sample in sample_lst:
    data_sample = data[data['sample'] == sample].copy().reset_index(drop=True)
    if len(data_sample) >= 5:
        sample_updated.append(sample)
        sample_n.append(len(data_sample))
        data_radial = pd.DataFrame()
        data_radial['0.05'] = [data_sample['radial_curve_normalized'][i][0] for i in range(len(data_sample))]
        data_radial['0.15'] = [data_sample['radial_curve_normalized'][i][1] for i in range(len(data_sample))]
        data_radial['0.25'] = [data_sample['radial_curve_normalized'][i][2] for i in range(len(data_sample))]
        data_radial['0.35'] = [data_sample['radial_curve_normalized'][i][3] for i in range(len(data_sample))]
        data_radial['0.45'] = [data_sample['radial_curve_normalized'][i][4] for i in range(len(data_sample))]
        data_radial['0.55'] = [data_sample['radial_curve_normalized'][i][5] for i in range(len(data_sample))]
        data_radial['0.65'] = [data_sample['radial_curve_normalized'][i][6] for i in range(len(data_sample))]
        data_radial['0.75'] = [data_sample['radial_curve_normalized'][i][7] for i in range(len(data_sample))]
        data_radial['0.85'] = [data_sample['radial_curve_normalized'][i][8] for i in range(len(data_sample))]
        data_radial['0.95'] = [data_sample['radial_curve_normalized'][i][9] for i in range(len(data_sample))]

        percentage_higher = [len(data_radial[data_radial['0.05'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.15'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.25'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.35'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.45'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.55'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.65'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.75'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.85'] > 1])/len(data_radial),
                             len(data_radial[data_radial['0.95'] > 1])/len(data_radial)]

        data_heatmap.loc[len(data_heatmap.index)] = percentage_higher
        data_heatmap_mean.loc[len(data_heatmap_mean.index)] = data_radial.mean()

print(sample_updated)
print(sample_n)

# hierarchical Clustering
# AHC (agglomerative, bottom-up)
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.

plt.figure(figsize=(12, 9))
clusters = shc.linkage(data_heatmap_mean, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/%s_ahc.pdf' % (save_folder, group))
plt.close()

nodes = R['ivl']
sample_nodes = [sample_updated[int(i)] for i in nodes]
print(sample_nodes)

data_heatmap_sort = pd.DataFrame(columns=data_heatmap_mean.columns)
for i in range(len(data_heatmap_mean)):
    data_heatmap_sort.loc[len(data_heatmap_sort.index)] = \
        data_heatmap_mean.iloc[int(dat.list_invert(nodes)[i])]

# heat map
plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap.pdf' % (save_folder, group))
plt.close()

plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap_sort, cbar=0, linewidths=2, vmax=1.2, vmin=0.8, square=True, cmap='coolwarm',
                  yticklabels=False)
plt.savefig('%s/%s_heatmap_sort.pdf' % (save_folder, group))
plt.close()

plt.subplots(figsize=(12, 9))
ax1 = sns.heatmap(data_heatmap_mean, cbar=0, linewidths=2, vmax=1.2, vmin=0.8, square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_mean.pdf' % (save_folder, group))
plt.close()

print("DONE!")