import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import math
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset6_dataset3-radial/" % master_folder
output_dir = "%sfigures3/" % master_folder

samples = ['40-0_5', '40-1_5', '40-2_5', '40-5_5', '40-10_5', '40-25_5', '40-50_5', '40-75_5', '40-100_5', '40-200_5', '40-300_5', '40-400_5', '40-500_5', '40-1000_5',
           '40-2000_5', '40-3000_5', '40-4000_5', '40-5000_5']
cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(samples))
line_color = [cmap(i) for i in x]
sns.set_palette(sns.color_palette(line_color))

hue_order = samples
hue_order1 = ['background', 'DNAFISH']

df_r_name = 'cen-r_40'
df_r = pd.read_csv(("%sdf_r_%s.txt" % (data_dir, df_r_name)), na_values=['.'], sep='\t')

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
for k in range(len(hue_order)):
    data_r = df_r[df_r['sample'] == hue_order[k]].copy().reset_index(drop=True)

    data_heatmap = pd.DataFrame(columns=x[:up])
    for s in hue_order1:
        data_sample = data_r[(data_r['sample_category'] == s) & (data_r['seg'] == 1)].copy().reset_index(drop=True)
        total_data_sample = len(data_sample)
        print(total_data_sample)
        data_radial = []
        for i in range(len(x[:up])):
            if x[i] == 0.025:
                n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)])
            else:
                n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))])
            data_radial.append(n_data_radial * 1.0/total_data_sample)
        data_heatmap.loc[len(data_heatmap.index)] = data_radial
    data_heatmap.index = hue_order1
    data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

    """plt.subplots(figsize=(12, len(hue_order1)))
    ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_before_normalization_%s.pdf' % (output_dir, df_r_name, hue_order[k]))
    plt.show()"""

    data_heatmap_normalized = pd.DataFrame(columns=x[:up])
    temp = []
    for i in data_heatmap.columns:
        temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
    data_heatmap_normalized.loc[0] = temp
    data_heatmap_normalized.index = [hue_order1[1]]
    data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
    df_seg_normalized.loc[len(df_seg_normalized.index)] = temp

    """plt.subplots(figsize=(12, len(hue_order)))
    ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(),
                      vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_DNAFISH_seg_n%s_cell%s_%s.pdf' % (output_dir, figure_name, n_dilation, len(data), hue_order[k]))
    plt.show()"""

df_seg_normalized.index = hue_order
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                  vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s_heatmap_DNAFISH_seg.pdf' % (output_dir, df_r_name))
plt.show()

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    plt.plot(x, df_seg_normalized.iloc[k].tolist(), color=line_color[k], label='%s' % hue_order[k])
plt.axhline(y=1, linestyle='--', color='r')
plt.xlabel('relative_r')
plt.ylabel('radial_curve')
plt.legend(loc=3)
plt.savefig('%s%s_rc-new_DNAFISH_seg.pdf' % (output_dir, df_r_name))
plt.show()