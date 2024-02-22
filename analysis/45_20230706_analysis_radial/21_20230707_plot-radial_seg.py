import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230706_analysis_radial/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample = 'Colo320DM_acidFISH_lamin_3d'
sample = 'Colo320HSR_acidFISH_lamin_3d'

df = pd.read_csv('%s%s/%s_radial_CN.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df_cell = pd.read_csv('%s%s/%s_n4_full.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
feature = ['radial_curve_nuclear', 'radial_curve_edge_nuclear',
           'radial_curve_laminB', 'radial_curve_edge_laminB']
for f in feature:
    df_cell[f] = [dat.str_to_float(df_cell[f][i]) for i in range(len(df_cell))]

name_order = ['laminB1 IF', 'MYC DNA FISH', 'hoechst stain']
hue_order = ['radial_curve_laminB', 'radial_curve_DNAFISH', 'radial_curve_nuclear']
hue_order2 = ['radial_curve_edge_laminB', 'radial_curve_edge_DNAFISH', 'radial_curve_edge_nuclear']
line_colors = [(255/255, 66/255, 161/255), (29/255, 177/255, 0/255), (0/255, 118/255, 186/255)]
sns.set_palette(sns.color_palette(line_colors))
hue_order1 = ['bg', 'DNAFISH']

# heatmap
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
up = 40

# seg
data_r = df
data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
data_r_sort['sample_category'] = ['DNAFISH'] * len(data_r_sort)
data_r_bg = data_r.copy()
data_r_bg['seg'] = [1] * len(data_r_bg)
data_r_bg['sample_category'] = ['bg'] * len(data_r_bg)
data_r_final = pd.concat([data_r_sort, data_r_bg], axis=0)

data_heatmap = pd.DataFrame(columns=x[:up])
for s in hue_order1:
    data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
    total_data_sample = len(data_sample)
    data_radial = []
    for i in range(len(x[:up])):
        if x[i] == 0.0125:
            n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.025)])
        else:
            n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.0125)) & (data_sample['relative_r'] <= (x[i]+0.0125))])
        data_radial.append(n_data_radial * 1.0/total_data_sample)
    data_heatmap.loc[len(data_heatmap.index)] = data_radial
data_heatmap.index = hue_order1
data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

temp = []
for i in data_heatmap.columns:
    temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])

vmax_lst = []
vmin_lst = []
x = np.arange(0.0125, 1, 0.025)
column_lst = [str('%.3f' % i) for i in x]
data_heatmap = pd.DataFrame(columns=column_lst)
for s in range(len(hue_order)):
    if hue_order[s] != 'radial_curve_DNAFISH':
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [df_cell[hue_order[s]][i][j] for i in range(len(df_cell))]
        data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
        vmax_lst.append(max(data_radial.mean()))
        vmin_lst.append(min(data_radial.mean()))
    else:
        data_heatmap.loc[len(data_heatmap.index)] = temp
        vmax_lst.append(max(temp))
        vmin_lst.append(min(temp))
data_heatmap.index = name_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=max(vmax_lst), vmin=min(vmin_lst), square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=2.3, vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial_curve
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df_cell)
    mean_curve = data_heatmap.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s, n=%s' % (name_order[k], number_nuclear))
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0, 2.3])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

# heatmap
x = np.arange(59.25, 0, -1.5)
x_label = 'absolute r'
up = 40

# seg
data_heatmap = pd.DataFrame(columns=x[:up])
for s in hue_order1:
    data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
    total_data_sample = len(data_sample)
    data_radial = []
    for i in range(len(x[:up])):
        if x[i] == 0.75:
            n_data_radial = len(data_sample[(data_sample['absolute_r'] >= 0) & (data_sample['absolute_r'] <= 1.5)])
        else:
            n_data_radial = len(data_sample[(data_sample['absolute_r'] > (x[i]-0.75)) & (data_sample['absolute_r'] <= (x[i]+0.75))])
        data_radial.append(n_data_radial * 1.0/total_data_sample)
    data_heatmap.loc[len(data_heatmap.index)] = data_radial
data_heatmap.index = hue_order1
data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

temp = []
for i in data_heatmap.columns:
    temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])

x = np.arange(59.25, 0, -1.5)
column_lst = [str('%.3f' % i) for i in x]
data_heatmap = pd.DataFrame(columns=column_lst)
for s in range(len(hue_order)):
    if hue_order2[s] != 'radial_curve_edge_DNAFISH':
        data_radial = pd.DataFrame()
        for j in range(len(column_lst)):
            data_radial[column_lst[j]] = [df_cell[hue_order2[s]][i][-(j+1)] for i in range(len(df_cell))]
        data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
    else:
        data_heatmap.loc[len(data_heatmap.index)] = temp
data_heatmap.index = name_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(),
                  vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=2.3,
                  vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial_curve
x = np.arange(0.75, 60, 1.5)
x_label = 'absolute r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df)
    mean_curve = data_heatmap.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(df_cell[df_cell['cn_group'] == hue_order[k]])))
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0, 2.3])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()