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

batch1 = '5uM_48hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg',
            'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

data_flt = pd.DataFrame()
for feature in features:
    data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
    data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)
data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)

data_hit1 = data_flt

batch2 = '5uM_24hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch2, batch2), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch2, batch2), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

features = ['n_filtered', 'n_neg', 'n_pos', 'fov_hoechst', 'pos_vs_neg',
            'per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M',
            'per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

data_flt = pd.DataFrame()
for feature in features:
    data_temp = data[(data['log2fc_%s' % feature] <= -1) | (data['log2fc_%s' % feature] >= 1)]
    data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)
data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)

data_hit2 = data_flt

# common
label_common = data_hit1.merge(data_hit2, how='inner', on=['label'])['label'].tolist()
if 'DMSO(DMSO)' in label_common:
    label_common.remove('DMSO(DMSO)')
# total
label_total = list(set(data_hit1['label'].tolist() + data_hit2['label'].tolist()))
if 'DMSO(DMSO)' in label_total:
    label_total.remove('DMSO(DMSO)')

batch2 = '5uM_24hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch2, batch2), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch2, batch2), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data2_common = data[data['label'].isin(label_common)].copy().reset_index(drop=True)
data2_common_sort = data2_common.sort_values(by='log2fc_fov_hoechst').reset_index(drop=True)
data2_total = data[data['label'].isin(label_total)].copy().reset_index(drop=True)
data2_total_sort = data2_total.sort_values(by='log2fc_fov_hoechst').reset_index(drop=True)

batch1 = '5uM_48hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch1, batch1), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data1_common = data[data['label'].isin(label_common)].copy().reset_index(drop=True)
data1_common['sort'] = data2_common['log2fc_fov_hoechst']
data1_common_sort = data1_common.sort_values(by='sort').reset_index(drop=True)
data1_total = data[data['label'].isin(label_total)].copy().reset_index(drop=True)
data1_total['sort'] = data2_total['log2fc_fov_hoechst']
data1_total_sort = data1_total.sort_values(by='sort').reset_index(drop=True)

features_survival = ['fov_hoechst', 'n_filtered', 'n_neg', 'n_pos', 'pos_vs_neg']
features_neg = ['per_neg_G1', 'per_neg_G1S', 'per_neg_S', 'per_neg_G2M']
features_pos = ['per_pos_G1', 'per_pos_G1S', 'per_pos_S', 'per_pos_G2M']

dfs = [data1_common_sort, data1_total_sort, data2_common_sort, data2_total_sort]
batches = [batch1, batch1, batch2, batch2]
names = ['common', 'total', 'common', 'total']

for i in range(4):
    data_flt = dfs[i]
    batch = batches[i]
    name = names[i]
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
    plt.savefig('%s/%s/%s_heatmap_%s.pdf' % (output_dir, batch, batch, name))
    plt.show()
