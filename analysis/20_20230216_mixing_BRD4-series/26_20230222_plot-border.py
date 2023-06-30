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
bg = 5500

hue_order = [sample1_neg, sample1_pos]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s_1.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df1 = df1[df1['circ_nuclear'] > 0.8].copy().reset_index(drop=True)
df1['observed_background'] = (df1['mean_int_DNAFISH']*df1['area_nuclear']-df1['mean_int_ecDNA']*df1['total_area_ecDNA'])/(df1['area_nuclear'] - df1['total_area_ecDNA'])
df1['detection_rate'] = (df1['mean_int_ecDNA']*df1['total_area_ecDNA'])/(df1['mean_int_DNAFISH']*df1['area_nuclear'])
df1['total_int_DNAFISH'] = df1['mean_int_DNAFISH']*df1['area_nuclear']
df1['mean_int_DNAFISH_corrected'] = [i-bg for i in df1['mean_int_DNAFISH'].tolist()]
df1['mean_int_DNAFISH_corrected'] = [0 if i<0 else i for i in df1['mean_int_DNAFISH_corrected'].tolist()]
df1['total_int_DNAFISH_corrected'] = df1['mean_int_DNAFISH_corrected'] * df1['area_nuclear']
df1['log2_total_int_DNAFISH_corrected'] = np.log2(df1['total_int_DNAFISH_corrected'])
df1['mean_int_ecDNA_corrected'] = df1['mean_int_ecDNA'] - bg

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

border = pd.DataFrame()
border['sample'] = df['sample'].tolist() + df['sample'].tolist() + df['sample'].tolist() + df['sample'].tolist() + df['sample'].tolist()
border['border'] = [-2] * len(df) + [-1] * len(df) + [0] * len(df) + [1] * len(df) + [2] * len(df)
border['mean_int'] = df['mean_int_border_m2'].tolist() + df['mean_int_border_m1'].tolist() + df['mean_int_border0'].tolist() + df['mean_int_border_p1'].tolist() + df['mean_int_border_p2'].tolist()
border['mean_int'] = border['mean_int']-5500

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 9))
sns.barplot(data=border, x='border', y='mean_int', hue='sample')
plt.savefig('%s%s_border.pdf' % (output_dir, figure_name))
plt.show()

