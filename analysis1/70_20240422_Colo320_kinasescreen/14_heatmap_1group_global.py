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

batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']

for batch in batches:
    data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
    data = pd.concat([data, data_cc], axis=1)
    data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

    search = 'Axl'
    search_lst = ['BI 2536', 'Centrinone-B', 'GSK461364', 'Rigosertib']
    search_lst = []
    data_flt = data[data['target'].str.contains(search)].copy().reset_index(drop=True)
    print(data_flt['treatment'].tolist())
    if len(search_lst) > 0:
        data_flt = data_flt[data_flt['treatment'].isin(search_lst)].copy().reset_index(drop=True)

    data_flt_na = data_flt.copy()
    features = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S', 'log2fc_per_G2M', 'log2fc_cytokinesis']
    for feature in features:
        data_flt_na.loc[(data_flt_na['rep'] == 2) & (data_flt_na['mean_n_filtered'] < 200), feature] = np.nan

    features = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis']

    data_hp = pd.DataFrame()
    data_hp_na = pd.DataFrame()
    for feature in features:
        data_hp[feature] = data_flt[feature]
        data_hp_na[feature] = data_flt_na[feature]
    data_hp.index = data_flt['label']
    data_hp_na.index = data_flt_na['label']

    fig, ax = plt.subplots(figsize=(7, len(data_flt) * 0.3 + 1))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(data_hp, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
    plt.savefig('%s/%s/%s_heatmap_global-new_%s.pdf' % (output_dir, batch, batch, search))
    plt.show()

    fig, ax = plt.subplots(figsize=(7, len(data_flt) * 0.3 + 1))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(data_hp_na, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
    plt.savefig('%s/%s/%s_heatmap_global-new_%s_na.pdf' % (output_dir, batch, batch, search))
    plt.show()


