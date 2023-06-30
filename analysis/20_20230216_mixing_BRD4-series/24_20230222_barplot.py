import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
from scipy import stats
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

d_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
n_bin = len(d_range)-1
all_df = []  # hue_order[0]: 0, 1, 2, 3, 4; hue_order[1]: 5, 6, 7, 8, 9
mean_df = []
err_lst = []
for i in range(len(hue_order)):
    data_sample = df[df['sample'] == hue_order[i]].copy().reset_index(drop=True)
    temp = []
    for j in range(len(d_range)-1):
        temp_sample = data_sample[(data_sample['total_int_DNAFISH_corrected'] >= d_range[j]*1E8) & (data_sample['total_int_DNAFISH_corrected'] < d_range[j+1]*1E8)]['mean_int_ecDNA_corrected'].tolist()
        all_df.append(temp_sample)
        temp.append(np.mean(temp_sample))
        err_lst.append(np.std(temp_sample))
    mean_df.append(temp)

for i in range(n_bin):
    print(stats.ttest_ind(all_df[i], all_df[i+n_bin], equal_var=False))

df_temp = pd.DataFrame(mean_df)
df_temp.columns = ['%.1f' % elem for elem in d_range][1:]
print(df_temp.columns)
df_temp.index = hue_order

dis = pd.DataFrame()
dis['sample'] = list(hue_order) * (len(d_range)-1)
dis['total_int_DNAFISH_corrected'] = [0.1] * len(hue_order) + [0.3] * len(hue_order) + [0.5] * len(hue_order) + [0.7] * len(hue_order) + [0.9] * len(hue_order)
dis['mean_int_ecDNA_corrected'] = df_temp['0.2'].tolist() + df_temp['0.4'].tolist() + df_temp['0.6'].tolist() + df_temp['0.8'].tolist() + df_temp['1.0'].tolist()
print(dis.head())

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 9))
ax = sns.barplot(data=dis, x='total_int_DNAFISH_corrected', y='mean_int_ecDNA_corrected', hue='sample')
"""x_coords = [p.get_x() + 0.5*p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
plt.errorbar(x=x_coords, y=y_coords, yerr=err_lst, fmt="none", c="k")"""
"""p_value = ['**', '***', '**', 'ns', 'ns']  # *** <0.01, ** <0.05, ns > 0.05
total_int_DNAFISH_corrected_bin = [0.1, 0.3, 0.5, 0.7, 0.9]
for i in range(n_bin):
    x1, x2 = -0.4+i, 0.4+i   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    y, h, col = dis[dis['total_int_DNAFISH_corrected'] == total_int_DNAFISH_corrected_bin[i]]['mean_int_ecDNA'].max() + 200, 2, 'k'
    print(y)
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, p_value[i], ha='center', va='bottom', color=col)"""
plt.savefig('%s%s_barplot_mean_int_ecDNA_corrected_group-by_total_int_DNAFISH_corrected.pdf' % (output_dir, figure_name))
plt.show()

print()

