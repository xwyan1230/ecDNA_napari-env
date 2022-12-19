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

sample = 'H3K27me2me3'
n = 8
df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_IF',
           'radial_curve_DNAFISH_normalized', 'radial_curve_IF_normalized',
           'nuclear_int', 'DNAFISH_int', 'IF_int', 'ecDNA_label',
           'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'

data = df
number_nuclear = len(data)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_DNAFISH'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_IF'].tolist())

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_DNAFISH'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_nuclear'][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_IF'][i], alpha=0.05, color=[line_colors[2][j]+0.05 for j in range(len(line_colors[2]))])
plt.plot(x, mean_curve1, color=line_colors[1], label='%s, n=%s' % ('ecDNA (MYC DNA FISH)', number_nuclear))
plt.plot(x, mean_curve2, color=line_colors[0], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x, mean_curve3, color=line_colors[2], label='%s, n=%s' % (sample, number_nuclear))
plt.plot(x, ci_lower1, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower3, color=line_colors[2], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[2], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()

mean_curve4, ci_lower4, ci_higher4 = dat.mean_list(data['radial_curve_DNAFISH_normalized'].tolist())
mean_curve5, ci_lower5, ci_higher5 = dat.mean_list(data['radial_curve_IF_normalized'].tolist())

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_DNAFISH_normalized'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_IF'][i], alpha=0.05, color=[line_colors[2][j] + 0.05 for j in range(len(line_colors[2]))])
plt.plot(x, mean_curve4, color=line_colors[1], label='%s, n=%s' % ('ecDNA (normalized)', number_nuclear))
plt.plot(x, ci_lower4, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher4, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, mean_curve3, color=line_colors[2], label='%s, n=%s' % (sample, number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[2], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[2], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()

# heatmap
data_heatmap = pd.DataFrame(columns=x)
hue_order = ['ecDNA', sample]
heatmap_order = ['radial_curve_DNAFISH_normalized', 'radial_curve_IF']
for s in heatmap_order:
    data_sample = df
    data_radial = pd.DataFrame()
    for i in range(len(x)):
        data_radial[x[i]] = [data_sample[s][j][i] for j in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap[:1], cbar=0, linewidths=2, vmax=data_heatmap[:1].values.max(), vmin=data_heatmap[:1].values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_ecDNA_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()