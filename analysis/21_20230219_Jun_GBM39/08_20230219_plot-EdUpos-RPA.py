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
sample1 = 'gbm39ec con'
sample2 = 'gbm39hsr con'
figure_name = 'gbm39ec_vs_hsr'
EdU_threshold = 30

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
df = df[df['mean_int_EdU'] > EdU_threshold].copy().reset_index()

# scatter plot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_RPA', x='mean_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_mean_int_RPA_vs_mean_int_DNAFISH_EdU%s.pdf' % (output_dir, figure_name, EdU_threshold))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='count_RPA', x='mean_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_count_RPA_vs_mean_int_DNAFISH_EdU%s.pdf' % (output_dir, figure_name, EdU_threshold))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='total_area_RPA', x='mean_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_total_area_RPA_vs_mean_int_DNAFISH_EdU%s.pdf' % (output_dir, figure_name, EdU_threshold))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
for i in range(len(hue_order)):
    df_temp = df[df['sample'] == hue_order[i]].copy().reset_index(drop=True)
    sns.scatterplot(data=df_temp, y='count_RPA', x='mean_int_DNAFISH', c=df_temp['mean_int_EdU'].tolist(), size=5)
    plt.legend()
    plt.savefig('%s/%s_count_RPA_vs_mean_int_DNAFISH_color_by_EdU_%s_EdU%s.pdf' % (output_dir, figure_name, hue_order[i], EdU_threshold))
    plt.show()
