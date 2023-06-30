import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
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
sample1 = 'DM-Ctrl_mix_mCh-Ctrl_forBRD3'
figure_name = sample1
pos_threshold = 20000
neg_threshold = 11000
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'

hue_order = [sample1_neg, sample1_pos]
DM_mean = 8445.679627934287
DM_mCh_mean = 8694.040717295931
hue_order_nor = [DM_mean/DM_mCh_mean if 'DM H2B-mCherry' in x else 1 for x in hue_order]
print(hue_order_nor)
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
df_sample = df[df['sample'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
df['normalized'] = [df[feature][i] * hue_order_nor[hue_order.index(df['sample'][i])] for i in range(len(df))]

sns.set_palette(sns.color_palette(line_colors))
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df, x='sample', y='normalized', order=hue_order)
sinaplot(data=df, x='sample', y='normalized', order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s_normalized_%s_n%s.pdf' % (output_dir, figure_name, feature, n_dilation))
plt.show()

for i in hue_order:
    print(i)
    temp = df[df['sample'] == i].copy().reset_index(drop=True)
    print(np.mean(temp['normalized']))