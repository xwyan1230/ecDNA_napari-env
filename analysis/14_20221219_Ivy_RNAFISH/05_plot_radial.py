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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'PC9'
n = 0
df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_RNAFISH', 'radial_curve_normalized']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20

data = df
number_nuclear = len(data)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_RNAFISH'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x[:up], data['radial_curve_RNAFISH'][i][:up], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
for i in range(len(data)):
    plt.plot(x[:up], data['radial_curve_nuclear'][i][:up], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
plt.plot(x[:up], mean_curve1[:up], color=line_colors[1], label='%s, n=%s' % ('nascent MYC RNA', number_nuclear))
plt.plot(x[:up], mean_curve2[:up], color=line_colors[0], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x[:up], ci_lower1[:up], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x[:up], ci_higher1[:up], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x[:up], ci_lower2[:up], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x[:up], ci_higher2[:up], color=line_colors[0], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x[:up], data['radial_curve_RNAFISH'][i][:up], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
plt.plot(x[:up], mean_curve1[:up], color=line_colors[1], label='%s, n=%s' % ('nascent MYC RNA', number_nuclear))
plt.plot(x[:up], ci_lower1[:up], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x[:up], ci_higher1[:up], color=line_colors[1], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.2, 2])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_RNAFISH_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()

# heatmap
data_heatmap = pd.DataFrame(columns=x[:up])
hue_order = ['nascent MYC RNA']
heatmap_order = ['radial_curve_RNAFISH']
for s in heatmap_order:
    data_sample = df
    data_radial = pd.DataFrame()
    for i in range(len(x[:up])):
        data_radial[x[i]] = [data_sample[s][j][i] for j in range(len(data_sample))]
    data_heatmap.loc[len(data_heatmap.index)] = data_radial.mean()
data_heatmap.index = hue_order

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_RNAFISH_%s_n%s.pdf' % (output_dir, sample, n))
plt.show()