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
sample1 = 'DM-JQ1_mix_mCh-Ctrl'
figure_name = sample1
pos_threshold = 10000
neg_threshold = 6000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM JQ1'

hue_order = [sample1_pos, sample1_neg]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
feature = 'mean_int_DNAFISH'

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
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

df = df1
df['r'] = np.sqrt(df['area_nuclear']/math.pi)
df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA']/math.pi)
df['total_area_ecDNA_sqrt_normalized'] = df['total_area_ecDNA_sqrt']/df['r']
df['total_area_ecDNA_normalized'] = df['total_area_ecDNA']/df['area_nuclear']
df_sample = df[df['sample'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample

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
