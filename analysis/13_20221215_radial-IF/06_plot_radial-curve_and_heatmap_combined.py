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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# sample_lst = ['laminB1', 'H3K9me3', 'H3K27me2me3', 'H3K9me2', 'H3K27Ac', 'H3K4me3']
sample_lst = ['laminB1', 'H3K9me3', 'H3K27me2me3', 'H3K27Ac', 'H3K4me3']
name_extension = 'all_but_H3K9me2'
n = 8
df_all = pd.DataFrame()
for sample in sample_lst:
    df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')
    df['sample'] = [sample] * len(df)
    df_all = pd.concat([df_all, df])

df_all = df_all.reset_index(drop=True)

cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/7)
line_color = [cmap(i) for i in x]
# color map
# laminB1: 0
# H3K9me3: 1
# H3K27me2me3: 2
# H3K9me2: 3
# H3K27Ac: 5
# H3K4me3: 6
line_colors = list(map(line_color.__getitem__, [0, 1, 2, 5, 6]))

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_IF',
           'radial_curve_DNAFISH_normalized', 'radial_curve_IF_normalized',
           'nuclear_int', 'DNAFISH_int', 'IF_int', 'ecDNA_label',
           'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df_all[f] = [dat.str_to_float(df_all[f][i]) for i in range(len(df_all))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x = ['%.3f' % elem for elem in x]
x_label = 'relative r'

plt.subplots(figsize=(12, 9))
for s in range(len(sample_lst)):
    data = df_all[df_all['sample'] == sample_lst[s]].copy().reset_index(drop=True)

    modified_color = [line_colors[s][j] + 0.1 for j in range(len(line_colors[s]))]
    modified_color = [i if i < 1 else 1 for i in modified_color]

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_IF'][i], alpha=0.05, color=modified_color)

for s in range(len(sample_lst)):
    data = df_all[df_all['sample'] == sample_lst[s]].copy().reset_index(drop=True)
    number_nuclear = len(data)

    mean_curve4, ci_lower4, ci_higher4 = dat.mean_list(data['radial_curve_DNAFISH_normalized'].tolist())
    mean_curve5, ci_lower5, ci_higher5 = dat.mean_list(data['radial_curve_IF'].tolist())
    modified_color = [line_colors[s][j] + 0.1 for j in range(len(line_colors[s]))]
    modified_color = [i if i < 1 else 1 for i in modified_color]
    plt.plot(x, mean_curve4, color=modified_color, linestyle='--')
    plt.plot(x, mean_curve5, color=line_colors[s], linewidth=3, label='%s, n=%s' % (sample_lst[s], number_nuclear))

plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_n%s_%s.pdf' % (output_dir, n, name_extension))
plt.show()

plt.subplots(figsize=(12, 9))
for s in range(len(sample_lst)):
    data = df_all[df_all['sample'] == sample_lst[s]].copy().reset_index(drop=True)

    modified_color = [line_colors[s][j] + 0.1 for j in range(len(line_colors[s]))]
    modified_color = [i if i < 1 else 1 for i in modified_color]

    for i in range(len(data)):
        plt.plot(x, data['radial_curve_DNAFISH_normalized'][i], alpha=0.02, color='blue')
    for i in range(len(data)):
        plt.plot(x, data['radial_curve_IF'][i], alpha=0.05, color=modified_color)

mean_curve4, ci_lower4, ci_higher4 = dat.mean_list(df_all['radial_curve_DNAFISH_normalized'].tolist())
plt.plot(x, mean_curve4, color='blue', linewidth=3, label='ecDNA, n=%s' % len(df_all))

for s in range(len(sample_lst)):
    data = df_all[df_all['sample'] == sample_lst[s]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve5, ci_lower5, ci_higher5 = dat.mean_list(data['radial_curve_IF'].tolist())
    plt.plot(x, mean_curve5, color=line_colors[s], linewidth=3, label='%s, n=%s' % (sample_lst[s], number_nuclear))

plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_ecDNA_combined_n%s_%s.pdf' % (output_dir, n, name_extension))
plt.show()

# heatmap
data_heatmap = pd.DataFrame(columns=x)
hue_order = ['ecDNA'] + sample_lst

for s in hue_order:
    if s == 'ecDNA':
        data_sample = df_all
        data_radial = pd.DataFrame()
        for i in range(len(x)):
            data_radial[x[i]] = [data_sample['radial_curve_DNAFISH_normalized'][j][i] for j in range(len(data_sample))]
    else:
        data_sample = df_all[df_all['sample'] == s].copy().reset_index(drop=True)
        data_radial = pd.DataFrame()
        for i in range(len(x)):
            data_radial[x[i]] = [data_sample['radial_curve_IF'][j][i] for j in range(len(data_sample))]

    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_n%s_%s.pdf' % (output_dir, n, name_extension))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap[:1], cbar=0, linewidths=2, vmax=data_heatmap[:1].values.max(), vmin=data_heatmap[:1].values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_ecDNA_n%s_%s.pdf' % (output_dir, n, name_extension))
plt.show()