import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
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
sample = 'C2'

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]

df = pd.read_csv('%s/%s/%s_n4_simplified.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.2].copy().reset_index(drop=True)
df_sample = df_sort[df_sort['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

print(len(df))

# heatmap
column_lst = ['0.025', '0.075', '0.125', '0.175', '0.225', '0.275', '0.325', '0.375', '0.425', '0.475', '0.525',
              '0.575', '0.625', '0.675', '0.725', '0.775', '0.825', '0.875', '0.925', '0.975']

data_heatmap = pd.DataFrame(columns=column_lst)

for s in hue_order:
    data_sample = df[df['group'] == s].copy().reset_index(drop=True)
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/%s/%s_heatmap_DNAFISH_n%s_overpoint2.pdf' % (output_dir, sample, sample, n_dilation))
plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.4, 1.6])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%s%s/%s_radial_curve_DNAFISH_n%s_overpoint2.pdf' % (output_dir, sample, sample, n_dilation))
plt.show()

