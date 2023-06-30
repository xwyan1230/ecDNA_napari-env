import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
from scipy import stats
import seaborn as sns
from shared.sinaplot import sinaplot
import numpy as np
import math
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
feature = 'total_int_DNAFISH_corrected'

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

for i in hue_order:
    print(i)
    temp = df[df['sample'] == i].copy().reset_index(drop=True)
    print(np.mean(temp[feature]))
t_stat = stats.ttest_ind(df[df['sample'] == hue_order[0]][feature].tolist(), df[df['sample'] == hue_order[1]][feature].tolist(), equal_var=False)
print(t_stat)
if t_stat[1] < 0.001:
    significance = '***'
elif t_stat[1] < 0.05:
    significance = '**'
else:
    significance = 'ns'

sns.set_palette(sns.color_palette(line_colors))
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df, x='sample', y=feature, order=hue_order)
sinaplot(data=df, x='sample', y=feature, order=hue_order, violin=False, scale='area')
x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df1[feature].max(), 2, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, significance, ha='center', va='bottom', color=col)
plt.savefig('%s/%s_%s_n%s.pdf' % (output_dir, figure_name, feature, n_dilation))
plt.show()

