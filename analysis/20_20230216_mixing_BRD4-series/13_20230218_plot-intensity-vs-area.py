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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

n_dilation = 4

# samples
sample1 = 'DM-Ctrl_mix_mCh-BRD4'
figure_name = sample1
pos_threshold = 21000
neg_threshold = 16000
sample1_pos = 'DM H2B-mCherry BRD4ko'
sample1_neg = 'DM'
new_seg = 11000
bg_neg = 4945.8102602298695
bg_pos = 5961.066231416957

hue_order = [sample1_neg, sample1_pos]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s_1.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df1 = df1[df1['circ_nuclear'] > 0.8].copy().reset_index(drop=True)

sample_lst = []
for i in range(len(df1)):
    if df1['mCherry_mean'].tolist()[i] < neg_threshold:
        sample_lst.append(sample1_neg)
    elif df1['mCherry_mean'].tolist()[i] > pos_threshold:
        sample_lst.append(sample1_pos)
    else:
        sample_lst.append('NA')
df1['sample'] = sample_lst
df1 = df1[df1['sample'].isin(hue_order)].copy().reset_index(drop=True)

df1['mean_int_DNAFISH_corrected'] = [df1['mean_int_DNAFISH'][i]-bg_neg if df1['sample'][i] == sample1_neg else df1['mean_int_DNAFISH'][i]-bg_pos for i in range(len(df1))]
df1['mean_int_DNAFISH_corrected'] = [0 if i<0 else i for i in df1['mean_int_DNAFISH_corrected'].tolist()]
df1['total_int_DNAFISH_corrected'] = df1['mean_int_DNAFISH_corrected'] * df1['area_nuclear']
df1['total_int_DNAFISH'] = df1['mean_int_DNAFISH'] * df1['area_nuclear']
df1['mean_int_ecDNA_corrected'] = [df1['mean_int_ecDNA'][i]-bg_neg if df1['sample'][i] == sample1_neg else df1['mean_int_ecDNA'][i]-bg_pos for i in range(len(df1))]
df1['mean_int_ecDNA_corrected'] = [0 if i<0 else i for i in df1['mean_int_ecDNA_corrected'].tolist()]
df1['observed_background'] = (df1['mean_int_DNAFISH']*df1['area_nuclear']-df1['mean_int_ecDNA']*df1['total_area_ecDNA'])/(df1['area_nuclear'] - df1['total_area_ecDNA'])
df1['observed_background_corrected'] = [df1['observed_background'][i]-bg_neg if df1['sample'][i] == sample1_neg else df1['observed_background'][i]-bg_pos for i in range(len(df1))]
df1['observed_background_corrected'] = [0 if i<0 else i for i in df1['observed_background_corrected'].tolist()]
df1['r'] = np.sqrt(df1['area_nuclear']/math.pi)
df1['total_area_ecDNA_sqrt'] = np.sqrt(df1['total_area_ecDNA']/math.pi)
df1['total_area_ecDNA_sqrt_normalized'] = df1['total_area_ecDNA_sqrt']/df1['r']
df1['total_area_ecDNA_normalized'] = df1['total_area_ecDNA']/df1['area_nuclear']
df = df1

# scatter plot
sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='total_area_ecDNA_normalized', x='mean_int_DNAFISH_corrected', hue='sample', hue_order=hue_order, s=5, alpha=0.8)
plt.savefig('%s/%s_total_area_ecDNA_normalized_vs_mean_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='total_area_ecDNA', x='total_int_DNAFISH_corrected', hue='sample', hue_order=hue_order, s=5, alpha=0.8)
plt.savefig('%s/%s_total_ecDNA_area_vs_total_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df[df['sample'] == hue_order[0]].copy().reset_index(drop=True), y='total_area_ecDNA_normalized', x='mean_int_DNAFISH_corrected', s=5, alpha=0.8)
plt.axhline(y=0.05**2, linestyle='--', color='r', label='0.05')
plt.legend()
plt.savefig('%s/%s_total_area_ecDNA_normalized_vs_mean_int_DNAFISH_corrected_control.pdf' % (output_dir, figure_name))
plt.show()

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='observed_background_corrected', x='total_int_DNAFISH_corrected', hue='sample', hue_order=hue_order, s=5, alpha=0.8)
plt.savefig('%s/%s_observed_background_corrected_vs_total_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_ecDNA_corrected', x='total_int_DNAFISH_corrected', hue='sample', hue_order=hue_order, s=3, alpha=0.8)
plt.savefig('%s/%s_mean_int_ecDNA_corrected_vs_total_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_ecDNA', x='total_int_DNAFISH', hue='sample', hue_order=hue_order, s=3, alpha=0.8)
plt.savefig('%s/%s_mean_int_ecDNA_vs_total_int_DNAFISH.pdf' % (output_dir, figure_name))
plt.show()

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.scatterplot(data=df, y='mean_int_ecDNA', x='total_int_DNAFISH_corrected', hue='sample', hue_order=hue_order, s=3, alpha=0.8)
plt.savefig('%s/%s_mean_int_ecDNA_vs_total_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()
