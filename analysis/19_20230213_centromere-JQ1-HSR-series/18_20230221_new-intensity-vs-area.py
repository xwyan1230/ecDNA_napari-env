import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import math
from shared.sinaplot import sinaplot
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230206_analysis_centromere-JQ1-DM-HSR/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 2

# samples
sample1 = 'DM-mCh_mix_HSR'
figure_name = sample1
pos_threshold = 10000
neg_threshold = 5500
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'HSR'

hue_order = [sample1_pos, sample1_neg]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df1 = df1[df1['circ_nuclear'] > 0.8].copy().reset_index(drop=True)
df1['observed_background'] = (df1['mean_int_DNAFISH']*df1['area_nuclear']-df1['mean_int_ecDNA']*df1['total_area_ecDNA'])/(df1['area_nuclear'] - df1['total_area_ecDNA'])
df1['detection_rate'] = (df1['mean_int_ecDNA']*df1['total_area_ecDNA'])/(df1['mean_int_DNAFISH']*df1['area_nuclear'])
df1['total_int_DNAFISH'] = df1['mean_int_DNAFISH']*df1['area_nuclear']

sample_lst = []
for i in range(len(df1)):
    if df1['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(sample1_neg)
    elif df1['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(sample1_pos)
    else:
        sample_lst.append('NA')
df1['sample'] = sample_lst

df = df1
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['total_area_ecDNA_normalized'] = df['total_area_ecDNA']/df['area_nuclear']
df_sample = df[df['sample'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

# scatter plot
"""sns.set_palette(sns.color_palette(line_colors))
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
plt.show()"""

"""sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
df_temp = df[df['sample'] == hue_order[1]].copy().reset_index(drop=True)
sns.scatterplot(data=df_temp, y='total_area_ecDNA_normalized', x='mean_int_DNAFISH', c=df_temp['area_nuclear'].tolist(), size=5, alpha=0.8)
plt.savefig('%s/%s_total_area_ecDNA_normalized_vs_mean_int_DNAFISH_color_by_area_nuclear.pdf' % (output_dir, figure_name))
plt.show()"""

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='observed_background', x='total_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_observed_background_vs_total_int_DNAFISH.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_ecDNA', x='total_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_mean_int_ecDNA_vs_total_int_DNAFISH.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_ecDNA', x='total_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.xlim([0, 4.5E8])
plt.savefig('%s/%s_mean_int_ecDNA_vs_total_int_DNAFISH_scale.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='total_area_ecDNA', x='total_int_DNAFISH', hue='sample', hue_order=hue_order, size=5, alpha=0.8)
plt.savefig('%s/%s_total_ecDNA_area_vs_total_int_DNAFISH.pdf' % (output_dir, figure_name))
plt.show()

