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

batch = '5uM_24hr'
n_feature_group = 2
cal_type = 'notinB'  # accepts "and" or "notinB"

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

"""features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg',
            'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']"""
features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg']
features1 = ['per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

data_flt = pd.DataFrame()
for feature in features:
    data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
    data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)

if n_feature_group == 2:
    data_flt1 = pd.DataFrame()
    for feature in features1:
        data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
        data_flt1 = pd.concat([data_flt1, data_temp], axis=0).drop_duplicates().reset_index(drop=True)
    if cal_type == 'and':
        data_flt = data_flt.merge(data_flt1, how='inner').reset_index(drop=True)
    else:
        data_flt = data_flt[~data_flt['label'].isin(data_flt1['label'].tolist())].copy().reset_index(drop=True)

data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)

features_survival = ['fov_hoechst', 'n_filtered', 'n_neg', 'n_pos', 'pos_vs_neg']
features_neg = ['per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M']
features_pos = ['per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

data_hp = pd.DataFrame()
for feature in features_survival:
    data_hp['log2fc_%s' % feature] = data_flt['log2fc_%s' % feature]
data_hp['rep_neg'] = data_flt['rep_neg']
data_hp['mean_n_neg_larger_than_100'] = data_flt['mean_n_neg_larger_than_100']
for feature in features_neg:
    data_hp['log2fc_%s' % feature] = data_flt['log2fc_%s' % feature]
data_hp['rep_pos'] = data_flt['rep_pos']
data_hp['mean_n_pos_larger_than_100'] = data_flt['mean_n_pos_larger_than_100']
for feature in features_pos:
    data_hp['log2fc_%s' % feature] = data_flt['log2fc_%s' % feature]
data_hp.index = data_flt['label']

fig, ax = plt.subplots(figsize=(7, len(data_flt) * 0.3 + 1))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_hit_survival_not_cellcycle.pdf' % (output_dir, batch, batch))
plt.show()


