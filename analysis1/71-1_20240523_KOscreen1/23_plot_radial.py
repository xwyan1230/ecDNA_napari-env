import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import shared.dataframe as dat
from scipy.stats import ttest_ind
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

sample = 'G2'
hueorder = ['mCherry', 'GFP']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
hue_order1 = ['bg', 'DNAFISH']

data = pd.read_csv('%s/%s/22_%s_radial_calculated.txt' % (output_dir, sample, sample), na_values=['.'], sep='\t')
df = data[data['group'].isin(['GFP', 'mCherry'])].copy().reset_index(drop=True)

# heatmap
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hueorder)):
    data_r = df[(df['group'] == hueorder[k])].copy().reset_index(drop=True)
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
df_seg_normalized.index = hueorder
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

df_seg_normalized.to_csv('%s%s/23_%s_radial_summary_relativer.txt' % (output_dir, sample, sample), index=False, sep='\t')

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s%s/23_%s_heatmap_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=1.8,
                  vmin=0.2, square=True, cmap='coolwarm')
plt.savefig('%s%s/23_%s_heatmap_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# radial_curve
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
plt.subplots(figsize=(12, 9))
for k in range(len(hueorder)):
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s' % hueorder[k])
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/23_%s_radial_curve_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.ylim([0.2, 1.8])
plt.savefig('%s%s/23_%s_radial_curve_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# heatmap
# x = np.arange(98.25, 0, -2.5)
x = np.arange(39.5, 0, -1)
x_label = 'absolute r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hueorder)):
    data_r = df[(df['group'] == hueorder[k])].copy().reset_index(drop=True)
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
df_seg_normalized.index = hueorder
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

df_seg_normalized.to_csv('%s%s/23_%s_radial_summary_absoluter.txt' % (output_dir, sample, sample), index=False, sep='\t')

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s%s/23_%s_heatmap_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.close()

plt.subplots(figsize=(12, len(hueorder)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=1.8,
                  vmin=0.2, square=True, cmap='coolwarm')
plt.savefig('%s%s/23_%s_heatmap_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# radial_curve
x = np.arange(1.25, 100, 2.5)
x_label = 'absolute r'
plt.subplots(figsize=(12, 9))
for k in range(len(hueorder)):
    mean_curve = df_seg_normalized.iloc[k]
    plt.plot(x, mean_curve, color=line_colors[k], label='%s' % hueorder[k])
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
# plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/23_%s_radial_curve_edge_DNAFISH_seg.pdf' % (output_dir, sample, sample))
plt.ylim([0.2, 1.8])
plt.savefig('%s%s/23_%s_radial_curve_edge_DNAFISH_seg_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

print("DONE!")