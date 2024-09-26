import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batch = 'point5uM_48hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
print(data.head())

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='log2fc_n_filtered', y='log2fc_fov_hoechst', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_n_filtered', y='log2fc_fov_hoechst', alpha=0.5, s=30, c='r')
data_flt = data[(data['log2fc_n_filtered'] < -2.5) & (data['log2fc_fov_hoechst'] > -1.5)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2fc_n_filtered'][i] + 0.03, y=data_flt['log2fc_fov_hoechst'][i], s=data_flt['label'][i],
             size=6, color=(0 / 255, 191 / 255, 255 / 255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_n-filtered_vs_fov-hoechst_log2fc.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='mean_n_filtered', y='mean_fov_hoechst', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='mean_n_filtered', y='mean_fov_hoechst', alpha=0.5, s=30, c='r')
plt.savefig('%s/%s/%s_n-filtered_vs_fov-hoechst_2rep.pdf' % (output_dir, batch, batch))
plt.show()
