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

sample = 'Colo320DM_acidFISH_lamin_3d'

df = pd.read_csv('%s%s/%s_radial_CN.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df_cell = pd.read_csv('%s%s/%s_n4_full.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')

line_colors = [(86/255, 193/255, 255/255), (142/255, 164/255, 199/255), (168/255, 147/255, 163/255),
               (195/255, 130/255, 130/255), (255/255, 100/255, 78/100)]
sns.set_palette(sns.color_palette(line_colors))
hue_order = ['<17', '[17, 17.5)', '[17.5, 18)', '[18, 18.5)', '>18.5']
hue_order1 = ['bg', 'DNAFISH']

# heatmap
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hue_order)):
    data_r = df[(df['sample'] == hue_order[k])].copy().reset_index(drop=True)
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
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.index = hue_order
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_DNAFISH_seg_CN.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=2.3,
                  vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_DNAFISH_seg_CN_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial_curve
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df)
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(df_cell[df_cell['cn_group'] == hue_order[k]])))
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0, 2.3])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH_seg_CN.pdf' % (output_dir, sample, sample))
plt.close()

# heatmap
x = np.arange(59.25, 0, -1.5)
x_label = 'absolute r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hue_order)):
    data_r = df[(df['sample'] == hue_order[k])].copy().reset_index(drop=True)
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
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.index = hue_order
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_edge_DNAFISH_seg_CN.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=2.3,
                  vmin=0, square=True, cmap='coolwarm')
plt.savefig('%s%s/%s_heatmap_edge_DNAFISH_seg_CN_fixscale.pdf' % (output_dir, sample, sample))
plt.close()

# radial_curve
x = np.arange(0.75, 60, 1.5)
x_label = 'absolute r'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    number_nuclear = len(df)
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s, n=%s' % (hue_order[k], len(df_cell[df_cell['cn_group'] == hue_order[k]])))
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylim([0, 2.3])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_edge_DNAFISH_seg_CN.pdf' % (output_dir, sample, sample))
plt.close()