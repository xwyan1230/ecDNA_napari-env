import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221219_analysis_Ivy_RNAFISH/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = 'manual'

sample = 'HCT116'
n = 0
df = pd.read_csv(("%s%s_radial_RNAseg_n%s_%s.txt" % (data_dir, sample, n, folder)), na_values=['.'], sep='\t')
print(len(df))

line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['RNAseg_label', 'int_r_to_edge', 'int_r_to_center', 'int_relative_r']

for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]


def img_to_pixel_int(mask: np.array, img: np.array):
    index = [i for i, e in enumerate(mask) if e != 0]
    out = list(map(img.__getitem__, index))
    return out


df_r = pd.DataFrame()
sample_lst = []
relative_r_lst = []
for i in range(len(df)):
    sample_lst = sample_lst + ['background'] * len(df['int_relative_r'][i]) + ['MYC RNAFISH'] * len([j for j in df['RNAseg_label'][i] if j > 0])
    relative_r_lst = relative_r_lst + df['int_relative_r'][i] + img_to_pixel_int(df['RNAseg_label'][i], df['int_relative_r'][i])
df_r['sample'] = sample_lst
df_r['relative_r'] = relative_r_lst
len_bg = len(df_r[df_r['sample'] == 'background'])
len_sample = len(df_r[df_r['sample'] == 'MYC RNAFISH'])
df_r['weights'] = [1.0/len_bg if i == 'background' else 1.0/len_sample for i in df_r['sample']]

sns.set_palette(sns.color_palette(line_colors))
hue_order = ['background', 'MYC RNAFISH']
fig, ax = plt.subplots(figsize=(9, 6))
fig.subplots_adjust(left=0.2)
ax = sns.histplot(data=df_r, x='relative_r', hue='sample', hue_order=hue_order, multiple='dodge', bins=20, weights=df_r['weights'])
# ['count', 'frequency', 'density', 'probability', 'proportion', 'percent']
plt.savefig('%s/histplot_relative_r_%s_cell%s_%s.pdf' % (output_dir, sample, len(df), folder))
plt.show()

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20
data_heatmap = pd.DataFrame(columns=x[:up])
hue_order = ['background', 'MYC RNAFISH']
for s in hue_order:
    data_sample = df_r[df_r['sample'] == s].copy().reset_index(drop=True)
    total_data_sample = len(data_sample)
    data_radial = []
    for i in range(len(x[:up])):
        if x[i] == 0.025:
            n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)])
        else:
            n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))])
        data_radial.append(n_data_radial * 1.0/total_data_sample)
    data_heatmap.loc[len(data_heatmap.index)] = data_radial
data_heatmap.index = hue_order
data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=2, vmax=data_heatmap.values.max(), vmin=data_heatmap.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_RNAFISH_RNAseg_before_normalization_%s_n%s_cell%s_%s.pdf' % (output_dir, sample, n, len(df), folder))
plt.show()

data_heatmap_normalized = pd.DataFrame(columns=x[:up])
temp = []
for i in data_heatmap.columns:
    temp.append(data_heatmap[i]['MYC RNAFISH']/data_heatmap[i]['background'])
data_heatmap_normalized.loc[0] = temp
data_heatmap_normalized.index = ['MYC RNAFISH']
data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]

plt.subplots(figsize=(12, len(hue_order)))
ax1 = sns.heatmap(data_heatmap_normalized, cbar=0, linewidths=2, vmax=data_heatmap_normalized.values.max(), vmin=data_heatmap_normalized.values.min(), square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_RNAFISH_RNAseg_%s_n%s_cell%s_%s.pdf' % (output_dir, sample, n, len(df), folder))
plt.show()