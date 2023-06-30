import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'DM_3_49pos'
figure_name = 'DM-MYC_woall'
df = pd.read_csv('%s%s_summary.txt' % (data_dir, sample), na_values=['.'], sep='\t')

hue_order = ['17.375', '17.625', '17.875', '18.125', '18.375', '18.625', '18.875', 'all']
cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(hue_order))
line_colors = [cmap(i) for i in x]
"""# color map
# laminB1: 0
# H3K9me3: 1
# H3K27me2me3: 2
# H3K9me2: 3
# H3K27Ac: 5
# H3K4me3: 6
line_colors = list(map(line_color.__getitem__, [0, 1, 2, 5, 6]))"""

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_DNAFISH_seg']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_sample = df[df['cellcycle'].isin(['G1', 'S', 'G2'])].copy().reset_index(drop=True)
df = df_sample
df['total_int_MYC'] = df['area_nuclear_IF'] * df['mean_int_MYC']
df['ln_total_int_MYC'] = np.log(df['total_int_MYC'])

print(len(df))

# heatmap
column_lst = ['0.025', '0.075', '0.125', '0.175', '0.225', '0.275', '0.325', '0.375', '0.425', '0.475', '0.525',
              '0.575', '0.625', '0.675', '0.725', '0.775', '0.825', '0.875', '0.925', '0.975']

data_heatmap = pd.DataFrame(columns=column_lst)

for s in hue_order:
    if s == 'all':
        data_sample = df
    else:
        cutoff = float(s)
        data_sample = df[(df['ln_total_int_MYC'] >= cutoff-0.125) & (df['ln_total_int_MYC'] < cutoff+0.125)].copy().reset_index(drop=True)
        print(len(data_sample))
    data_radial = pd.DataFrame()
    for j in range(len(column_lst)):
        data_radial[column_lst[j]] = [data_sample['radial_curve_DNAFISH'][i][j] for i in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()

data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'

plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    if hue_order[k] == 'all':
        data = df
    else:
        cutoff = float(hue_order[k])
        data = df[(df['ln_total_int_MYC'] >= cutoff - 0.125) & (df['ln_total_int_MYC'] < cutoff + 0.125)].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_DNAFISH_n%s_%s.pdf' % (output_dir, n_dilation, figure_name))
plt.show()

