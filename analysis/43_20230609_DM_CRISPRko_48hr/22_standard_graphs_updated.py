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
sample = 'C3'

GFP_sample = 'GFP'
mCherry_sample = 'mCherry'
hue_order = [GFP_sample, mCherry_sample]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
sns.set_palette(sns.color_palette(line_colors))

df = pd.read_csv('%s/%s/%s_n4_updated1.txt' % (data_dir, sample, sample), na_values=['.'], sep='\t')
df['total_int_DNAFISH'] = df['mean_int_DNAFISH'] * df['area_nuclear']
feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'n_ecDNA', 'total_area_ratio_ecDNA', 'total_area_ecDNA']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df_sample = df[df['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
df_GFP = df[df['group'] == 'GFP'].copy().reset_index(drop=True)
df_mCherry = df[df['group'] == 'mCherry'].copy().reset_index(drop=True)

# n_ecDNA curve
"""x = np.arange(0, 65*1000, 1000)
x_label = 'intensity'
feature = 'n_ecDNA'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data[feature].tolist())
    for i in range(len(data)):
        plt.plot(x, data[feature][i], alpha=0.001, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 20])
plt.ylabel(feature)
plt.legend()
# plt.savefig('%s%s/%s_%s_updated.pdf' % (output_dir, sample, sample, feature))
plt.show()"""

x = np.arange(2500, 15*5000+2500, 5000)
x_label = 'intensity'
feature = 'total_area_ratio_ecDNA'
plt.subplots(figsize=(12, 9))
for k in range(len(hue_order)):
    data = df[df['group'] == hue_order[k]].copy().reset_index(drop=True)
    number_nuclear = len(data)
    mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data[feature].tolist())
    print(mean_curve3)
    print(sum(mean_curve3))
    for i in range(len(data)):
        plt.plot(x, data[feature][i], alpha=0.01, color=line_colors[k])
    plt.plot(x, mean_curve3, color=line_colors[k], label='%s, n=%s' % (hue_order[k], number_nuclear))
    plt.plot(x, ci_lower3, color=line_colors[k], linestyle='--', linewidth=0.5)
    plt.plot(x, ci_higher3, color=line_colors[k], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0, 1])
plt.ylabel(feature)
plt.legend()
# plt.savefig('%s%s/%s_%s_updated.pdf' % (output_dir, sample, sample, feature))
plt.show()
# plt.close()