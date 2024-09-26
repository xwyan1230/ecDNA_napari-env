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

# target = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % data_dir, na_values=['.'])
# print(target.head())

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data['divide'] = data['mean_fov_hoechst']/data['mean_n_filtered']
data['cc_ratio'] = (data['cc_score_pos']+1)/(data['cc_score_neg']+1)
mean_divide = np.mean(data['divide'])
data_sort = data[(data['mean_n_filtered'] >= 200) & (data['divide'] <= 1.6*mean_divide) & (data['divide'] >= 0.4*mean_divide)].copy().reset_index(drop=True)
data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
print(data.head())


def func(x, a):
    return a * x


plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data, x='mean_n_filtered', y='mean_fov_hoechst', alpha=0.5, s=30)
sns.scatterplot(data=data_sort, x='mean_n_filtered', y='mean_fov_hoechst', color='r', alpha=0.5, s=30)
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_n-neg_vs_n-pos_filter.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='mean_n_neg', y='mean_n_pos', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='mean_n_neg', y='mean_n_pos', color='r', alpha=0.5, s=30)
plt.xlim([-50, 2000])
plt.ylim([-50, 2000])
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_n-neg_vs_n-pos_sort.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='mean_n_neg', y='mean_n_pos', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='mean_n_neg', y='mean_n_pos', color='r', alpha=0.5, s=30)
data_flt = data_sort[(data_sort['log2fc_pos_vs_neg'] < -0.7) | (data_sort['log2fc_pos_vs_neg'] > 0.7)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['mean_n_neg'][i]+0.005, y=data_flt['mean_n_pos'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
plt.xlim([-50, 2000])
plt.ylim([-50, 2000])
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_n-neg_vs_n-pos_sort_label.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='cc_score_neg', y='cc_score_pos', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='cc_score_neg', y='cc_score_pos', color='r', alpha=0.5, s=30)
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_cc-neg_vs_cc-pos_sort.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='cc_score_neg', y='cc_score_pos', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='cc_score_neg', y='cc_score_pos', color='r', alpha=0.5, s=30)
data_flt = data_sort[((data_sort['cc_score_neg'] > 1) & (data_sort['cc_score_pos'] > 1)) & ((data_sort['cc_ratio'] < 0.4) | (data_sort['cc_ratio'] > 1.6))].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['cc_score_neg'][i]+0.005, y=data_flt['cc_score_pos'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_cc-neg_vs_cc-pos_sort_label.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_ratio', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_ratio', color='r', alpha=0.5, s=30)
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-ratio_sort.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_ratio', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_ratio', color='r', alpha=0.5, s=30)
data_flt = data_sort[(data_sort['log2fc_pos_vs_neg'] < -0.7) | (data_sort['log2fc_pos_vs_neg'] > 0.7)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2fc_pos_vs_neg'][i]+0.005, y=data_flt['cc_ratio'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-ratio_sort_label.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_score', color='r', alpha=0.5, s=30)
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-score_sort.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_score', color='r', alpha=0.5, s=30)
data_flt = data_sort[(data_sort['log2fc_pos_vs_neg'] < -0.7) | (data_sort['log2fc_pos_vs_neg'] > 0.7)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2fc_pos_vs_neg'][i]+0.005, y=data_flt['cc_score'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-score_sort_label.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_score', color='r', alpha=0.5, s=30)
plt.legend()
plt.ylim([-1, 10])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-score_sort_part.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_pos_vs_neg', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_pos_vs_neg', y='cc_score', color='r', alpha=0.5, s=30)
plt.ylim([-1, 10])
data_flt = data_sort[(data_sort['log2fc_pos_vs_neg'] < -0.7) | (data_sort['log2fc_pos_vs_neg'] > 0.7)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2fc_pos_vs_neg'][i]+0.005, y=data_flt['cc_score'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_survival-ratio_vs_cc-score_sort_label_part.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_fov_hoechst', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_fov_hoechst', y='cc_score', color='r', alpha=0.5, s=30)
plt.legend()
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_log2fc_fov_hoechst_vs_cc-score_sort.pdf' % (output_dir, batch, batch))
plt.show()

plt.subplots(figsize=(9, 9))
sns.scatterplot(data=data_sort, x='log2fc_fov_hoechst', y='cc_score', alpha=0.5, s=30)
sns.scatterplot(data=data_ctrl, x='log2fc_fov_hoechst', y='cc_score', color='r', alpha=0.5, s=30)
data_flt = data_sort[(data_sort['log2fc_fov_hoechst'] < -1)].copy().reset_index(drop=True)
for i in range(len(data_flt)):
    plt.text(x=data_flt['log2fc_fov_hoechst'][i]+0.005, y=data_flt['cc_score'][i]+0.005, s=data_flt['label'][i], size=7, color=(0/255, 191/255, 255/255))
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_log2fc_fov_hoechst_vs_cc-score_sort_label_part.pdf' % (output_dir, batch, batch))
plt.show()








