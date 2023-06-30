import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import os
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230504_analysis_DM_DNAFISH_ctrl/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

# new approach
hue_order1 = ['background', 'DNAFISH']
row_lst = ['B', 'C', 'D', 'E', 'F', 'G']
column_lst = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
sample_lst_total = []
for row in row_lst:
    sample_lst_total = sample_lst_total + ['%s%s' % (row, column) for column in column_lst]
ctrl_lst = ['C2', 'C10', 'D6', 'F2', 'F10']

# heatmap
x = np.arange(0.025, 1, 0.05)
x_label = 'relative r'
up = 20

# seg
for current_sample in sample_lst_total:
    df_seg_normalized = pd.DataFrame(columns=x[:up])
    sample_lst = ctrl_lst + [current_sample]
    for sample in sample_lst:
        if not os.path.exists('%s/txt_radial/%s_radial.txt' % (data_dir1, sample)):
            df_seg_normalized.loc[len(df_seg_normalized.index)] = [0] * up
        else:
            data_r = pd.read_csv('%s/txt_radial/%s_radial.txt' % (data_dir1, sample), na_values=['.'], sep='\t')
            data_heatmap = pd.DataFrame(columns=x[:up])
            for s in hue_order1:
                data_sample = data_r[(data_r['sample_category'] == s) & (data_r['seg'] == 1)].copy().reset_index(drop=True)
                total_data_sample = len(data_sample)
                data_radial = []
                for i in range(len(x[:up])):
                    if x[i] == 0.025:
                        n_data_radial = len(data_sample[(data_sample['relative_r'] >= 0) & (data_sample['relative_r'] <= 0.05)])
                    else:
                        n_data_radial = len(data_sample[(data_sample['relative_r'] > (x[i]-0.025)) & (data_sample['relative_r'] <= (x[i]+0.025))])
                    data_radial.append(n_data_radial * 1.0/total_data_sample)
                data_heatmap.loc[len(data_heatmap.index)] = data_radial
            data_heatmap.index = hue_order1
            data_heatmap.columns = ['%.3f' % elem for elem in x[:up]]

            data_heatmap_normalized = pd.DataFrame(columns=x[:up])
            temp = []
            for i in data_heatmap.columns:
                temp.append(data_heatmap[i][hue_order1[1]] / data_heatmap[i][hue_order1[0]])
            data_heatmap_normalized.loc[0] = temp
            data_heatmap_normalized.index = [hue_order1[1]]
            data_heatmap_normalized.columns = ['%.3f' % elem for elem in x[:up]]
            df_seg_normalized.loc[len(df_seg_normalized.index)] = temp

    df_seg_normalized.index = sample_lst
    df_seg_normalized.columns = ['%.3f' % elem for elem in x[:up]]

    plt.subplots(figsize=(12, len(sample_lst)))
    ax1 = sns.heatmap(df_seg_normalized, cbar=0, linewidths=2, vmax=df_seg_normalized.values.max(),
                      vmin=df_seg_normalized.values.min(), square=True, cmap='coolwarm')
    if not os.path.exists("%s/radial/" % output_dir):
        os.makedirs("%s/radial/" % output_dir)
    plt.savefig('%s/radial/%s_radial_seg.pdf' % (output_dir, current_sample))
    plt.close()

