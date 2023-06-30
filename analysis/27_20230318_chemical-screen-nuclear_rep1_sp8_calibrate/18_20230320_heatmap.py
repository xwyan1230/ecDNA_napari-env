import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
from scipy.stats import ks_2samp
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230318_analysis_chemical-screen-nuclear_rep1_sp8_calibrate/"
data_dir = "%sdata/" % master_folder
data_dir1 = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_nuclear_convex_dilation = -3

# samples
exp = 'DM_2hr'
ctrl = pd.read_csv('%s%s_ctrl.txt' % (data_dir, exp), na_values=['.'], sep='\t')
"""features = ['area_nuclear_mean_log2FC', 'mean_int_nuclear_cal_mean_log2FC', 'mean_int_MYC_cal_mean_log2FC',
            'total_int_nuclear_cal_mean_log2FC', 'total_int_MYC_cal_mean_log2FC']"""
features = ['area_nuclear_mean_log2FC', 'mean_int_nuclear_mean_log2FC', 'mean_int_MYC_mean_log2FC',
            'total_int_nuclear_mean_log2FC', 'total_int_MYC_mean_log2FC']
vmax_lst = [1.5, 1.5, 1.5, 1.5, 1.5]
vmin_lst = [-1.5, -1.5, -1.5, -1.5, -1.5]
ctrl['well'] = [ctrl['sample'][i].split('_')[0] for i in range(len(ctrl))]

# load data
df = pd.read_csv('%s%s_n%s_sum_1st.txt' % (data_dir1, exp, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
df_ref = pd.read_csv('%sDM_6hr_n%s_sum_1st.txt' % (data_dir1, n_nuclear_convex_dilation), na_values=['.'], sep='\t')
df['sort'] = df_ref['total_int_MYC_mean_log2FC']
df = df.sort_values(by='sort').reset_index(drop=True)

# scatter plot
for f in range(len(features)):
    data_frame = pd.DataFrame(columns=df['compound'].tolist())
    data_frame.loc[0] = df[features[f]].tolist()
    data_frame.index = ['%s_%s' % (exp, features[f])]
    plt.subplots(figsize=(len(df), 3))
    ax1 = sns.heatmap(data_frame, cbar=0, linewidths=2, vmax=vmax_lst[f], vmin=vmin_lst[f], square=True, cmap='coolwarm')
    plt.savefig('%s/%s_heatmap_%s.pdf' % (output_dir, exp, features[f]))
    plt.show()