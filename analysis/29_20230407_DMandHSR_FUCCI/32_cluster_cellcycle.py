import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import math
import matplotlib as mpl
import numpy as np
from shared.sinaplot import sinaplot
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230407_analysis_DMandHSR_FUCCI/"
data_dir = "%sdata/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'DM_3_49pos'
figure_name = 'DM-cellcycle'
df = pd.read_csv('%s%s_summary.txt' % (data_dir, sample), na_values=['.'], sep='\t')

hue_order = ['17.375', '17.625', '17.875', '18.125', '18.375', '18.625', '18.875']
hue_order1 = ['background', 'DNAFISH']
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(hue_order))
line_colors1 = [cmap(i) for i in x]

df_sample = df[df['cellcycle'].isin(['G1', 'S', 'G2'])].copy().reset_index(drop=True)
df = df_sample
df['total_int_MYC'] = df['area_nuclear_IF'] * df['mean_int_MYC']
df['ln_total_int_MYC'] = np.log(df['total_int_MYC'])
print(len(df))

df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['dis_to_hub_area_normalized'] = df['dis_to_hub_area']/df['r']

# data filter
# df_sort = df
df_sort = df[df['total_area_ecDNA_sqrt_normalized'] > 0.15].copy().reset_index(drop=True)
# df_sort = df_sample[(df_sample['total_area_ecDNA_sqrt_normalized'] > 0.1) & (df_sample['mean_int_DNAFISH'] <= 8500)].copy().reset_index(drop=True)
print(len(df_sort))

# scatter plot
"""sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA_sqrt_normalized', y='dis_to_hub_area_normalized', hue='sample', hue_order=hue_order)
plt.savefig('%s/%s_dis_to_hub_area_normalized_vs_total_area_sqrt_normalized.pdf' % (output_dir, figure_name))
plt.show()"""
hue_order2 = ['G1', 'S', 'G2']
fig, ax = plt.subplots(figsize=(2*len(hue_order2)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='MYC_group1', y='dis_to_hub_area_normalized', order=hue_order2)
sinaplot(data=df_sort, x='cellcycle', y='dis_to_hub_area_normalized', order=hue_order2, violin=False, scale='area')
plt.savefig('%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, figure_name))
plt.show()

"""fig, ax = plt.subplots(figsize=(2*len(hue_order2)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df_sort, x='sample', y='dis_to_hub_area_normalized', order=hue_order)
sns.barplot(data=df_sort, x='cellcycle', y='dis_to_hub_area_normalized', order=hue_order2)
# plt.savefig('%s/%s_dis_to_hub_area_normalized.pdf' % (output_dir, figure_name))
plt.show()"""