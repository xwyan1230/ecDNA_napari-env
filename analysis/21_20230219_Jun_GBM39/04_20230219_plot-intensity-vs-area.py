import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
import math

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230219_analysis_Jun_EGFR_RPAs33p_Edu/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 0

# samples
sample1 = 'gbm39ec hu'
sample2 = 'gbm39hsr hu'
figure_name = 'gbm39ec_vs_hsr_hu'

hue_order = [sample1, sample2]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df2 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample2, n_dilation)), na_values=['.'], sep='\t')
df1['sample'] = [sample1] * len(df1)
df2['sample'] = [sample2] * len(df2)
print(len(df1))
df1 = df1[df1['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
print(len(df1))
print(len(df2))
df2 = df2[df2['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
print(len(df2))

feature = ['percentage_area_curve_ecDNA', 'cum_area_ind_ecDNA_filled', 'g']
for f in feature:
    df1[f] = [dat.str_to_float(df1[f][i]) for i in range(len(df1))]
    df2[f] = [dat.str_to_float(df2[f][i]) for i in range(len(df2))]

df = pd.concat([df1, df2], axis=0)
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']
df['total_area_ecDNA_normalized'] = df['total_area_ecDNA']/df['area_nuclear']

# scatter plot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_ecDNA_normalized', y='mean_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_mean_int_DNAFISH_vs_total_area_ecDNA_normalized.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='total_area_ecDNA_normalized', x='mean_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_total_area_ecDNA_normalized_vs_mean_int_DNAFISH.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df[df['sample'] == hue_order[0]].copy().reset_index(drop=True), y='total_area_ecDNA_normalized', x='mean_int_DNAFISH', size=5, alpha=0.8)
plt.axhline(y=0.05**2, linestyle='--', color='r', label='0.05')
plt.legend()
plt.savefig('%s/%s_total_area_ecDNA_normalized_vs_mean_int_DNAFISH_control.pdf' % (output_dir, figure_name))
plt.show()