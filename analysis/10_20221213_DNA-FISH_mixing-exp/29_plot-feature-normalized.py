import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample1 = 'H2B+POLR3D'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 12000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM POLR3Dko'

hue_order = [sample1_pos, sample1_neg]
DM_mean = 13613.241793656312
DM_mCh_mean = 14217.03780483494
hue_order_nor = [DM_mean/DM_mCh_mean if 'DM H2B-mCherry' in x else 1 for x in hue_order]
print(hue_order_nor)
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
feature = 'mean_int_DNAFISH'

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')

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
df_sample = df[df['sample'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
df['normalized'] = [df[feature][i] * hue_order_nor[hue_order.index(df['sample'][i])] for i in range(len(df))]

sns.set_palette(sns.color_palette(line_colors))
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
sns.violinplot(data=df, x='sample', y='normalized', order=hue_order)
plt.savefig('%s/%s_normalized_%s_n%s.pdf' % (output_dir, figure_name, feature, n_dilation))
plt.show()