import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
from shared.sinaplot import sinaplot
import seaborn as sns
import numpy as np
import math
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'G7'

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
sns.set_palette(sns.color_palette(line_colors))

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df['total_int_DNAFISH'] = df['mean_int_DNAFISH'] * df['area_nuclear']
feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_edge_nuclear', 'radial_curve_edge_DNAFISH', 'percentage_area_curve_ecDNA', 'n_ecDNA_lst',
           'total_area_DNAFISH_lst', 'total_area_ratio_DNAFISH_lst']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df_sample = df[df['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

# radial heatmap
print("Plotting total radial curve...")
temp = np.arange(98.75, 0, -2.5)
column_lst = [str(i) for i in temp]
data_heatmap = pd.DataFrame(columns=column_lst)
for s in hue_order:
    data_sample = df[df['group'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_edge_DNAFISH'][i][-(j+1)] for i in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_edge_DNAFISH.pdf' % (output_dir, sample, sample))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=1.3, vmin=0.6, square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_edge_DNAFISH_fixscale.pdf' % (output_dir, sample, sample))
plt.show()

# radial curve
x = np.arange(98.75, 0, -2.5)
x_label = 'r to edge (reverse)'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_edge_DNAFISH'].tolist())
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_edge_DNAFISH'][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve_edge')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_edge_DNAFISH.pdf' % (output_dir, sample, sample))
plt.show()