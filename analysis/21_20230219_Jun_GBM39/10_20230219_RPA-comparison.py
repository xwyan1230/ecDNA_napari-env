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
sample2 = 'gbm39ec hu'
sample3 = 'gbm39hsr con'
sample4 = 'gbm39hsr hu'
EdU_threshold = 30

hue_order = [sample1, sample2, sample3, sample4]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25), 'b', 'g']  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df2 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample2, n_dilation)), na_values=['.'], sep='\t')
df3 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample3, n_dilation)), na_values=['.'], sep='\t')
df4 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample4, n_dilation)), na_values=['.'], sep='\t')
df1['sample'] = [sample1] * len(df1)
df2['sample'] = [sample2] * len(df2)
df3['sample'] = [sample3] * len(df3)
df4['sample'] = [sample4] * len(df4)
df1 = df1[df1['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
df2 = df2[df2['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
df3 = df3[df3['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
df4 = df4[df4['circ_nuclear'] > 0.7].copy().reset_index(drop=True)

df = pd.concat([df1, df2, df3, df4], axis=0).reset_index()
df = df[df['count_RPA']>0].copy().reset_index(drop=True)

zero_per = []
nonzero_per = []
for i in hue_order:
    df_temp = df[df['sample'] == i].copy().reset_index(drop=True)
    zero_per.append(len(df_temp[df_temp['count_RPA'] == 0])/len(df_temp))
    nonzero_per.append(1-len(df_temp[df_temp['count_RPA'] == 0])/len(df_temp))

# df = df[df['mean_int_EdU'] > EdU_threshold].copy().reset_index()

# barplot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='count_RPA', y='mean_int_RPA', hue='sample', hue_order=hue_order, size=5, alpha=0.7)
plt.savefig('%sRPA_count_vs_mean_int_RPA_nonzero.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='total_area_RPA', y='mean_int_RPA', hue='sample', hue_order=hue_order, size=5, alpha=0.7)
plt.savefig('%stotal_area_RPA_vs_mean_int_RPA_nonzero.pdf' % output_dir)
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, x='count_RPA', y='total_area_RPA', hue='sample', hue_order=hue_order, size=5, alpha=0.7)
plt.savefig('%sRPA_count_vs_total_area_RPA_nonzero.pdf' % output_dir)
plt.show()

