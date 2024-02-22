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
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_radial-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'laminB1'
n = 4

df = pd.read_csv(("%s%s_radial_r.txt" % (data_dir, sample)), na_values=['.'], sep='\t')
df_cell = pd.read_csv(("%s%s_radial_new_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

hue_order1 = ['bg', 'DNAFISH']

# heatmap
x = np.arange(0.0125, 1, 0.025)
x_label = 'relative r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
data_r = df
data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
data_r_sort['sample_category'] = ['DNAFISH'] * len(data_r_sort)
data_r_bg = data_r.copy()
data_r_bg['seg'] = [1] * len(data_r_bg)
data_r_bg['sample_category'] = ['bg'] * len(data_r_bg)
data_r_final = pd.concat([data_r_sort, data_r_bg], axis=0)

data_heatmap = pd.DataFrame(columns=x[:up])
for s in hue_order1:
    data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
    total_data_sample = len(data_sample)
    data_radial = []
    for i in range(len(x[:up])):
        if x[i] == 0.0125:
            n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.025)])
        else:
            n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.0125)) & (data_sample['relative_r'] <= (x[i]+0.0125))])
        data_radial.append(n_data_radial * 1.0/total_data_sample)
    data_heatmap.loc[len(data_heatmap.index)] = data_radial
data_heatmap.index = hue_order1
data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

temp = []
for i in data_heatmap.columns:
    temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]
df_seg_normalized.to_csv('%s/%s_radial_summary_relativer.txt' % (output_dir, sample), index=False, sep='\t')

# heatmap
x = np.arange(98.75, 0, -2.5)
x_label = 'absolute r'
up = 40
df_seg_normalized = pd.DataFrame(columns=x[:up])

# seg
data_r = df
data_r_sort = data_r[data_r['seg'] == 1].copy().reset_index(drop=True)
data_r_sort['sample_category'] = ['DNAFISH'] * len(data_r_sort)
data_r_bg = data_r.copy()
data_r_bg['seg'] = [1] * len(data_r_bg)
data_r_bg['sample_category'] = ['bg'] * len(data_r_bg)
data_r_final = pd.concat([data_r_sort, data_r_bg], axis=0)

data_heatmap = pd.DataFrame(columns=x[:up])
for s in hue_order1:
    data_sample = data_r_final[(data_r_final['sample_category'] == s) & (data_r_final['seg'] == 1)].copy().reset_index(drop=True)
    total_data_sample = len(data_sample)
    data_radial = []
    for i in range(len(x[:up])):
        if x[i] == 0.75:
            n_data_radial = len(data_sample[(data_sample['absolute_r'] >= 0) & (data_sample['absolute_r'] <= 1.5)])
        else:
            n_data_radial = len(data_sample[(data_sample['absolute_r'] > (x[i]-0.75)) & (data_sample['absolute_r'] <= (x[i]+0.75))])
        data_radial.append(n_data_radial * 1.0/total_data_sample)
    data_heatmap.loc[len(data_heatmap.index)] = data_radial
data_heatmap.index = hue_order1
data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

temp = []
for i in data_heatmap.columns:
    temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
df_seg_normalized.loc[len(df_seg_normalized.index)] = temp
df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

df_seg_normalized.to_csv('%s/%s_radial_summary_absoluter.txt' % (output_dir, sample), index=False, sep='\t')

print("DONE!")