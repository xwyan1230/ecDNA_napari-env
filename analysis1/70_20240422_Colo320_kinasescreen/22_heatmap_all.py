import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import os

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batch = '5uM_48hr'

features_survival = ['label', 'mean_n_filtered', 'log2fc_fov_hoechst', 'log2fc_n_filtered', 'log2fc_n_neg', 'log2fc_n_pos', 'log2fc_pos_vs_neg']
features_cellcycle = ['label', 'mean_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis', 'cc_score_neg', 'log2fc_per_neg_G1', 'log2fc_per_neg_G1S',
                           'log2fc_per_neg_S', 'log2fc_per_neg_G2M', 'cc_score_pos',
                           'log2fc_per_pos_G1', 'log2fc_per_pos_G1S', 'log2fc_per_pos_S', 'log2fc_per_pos_G2M']

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data_survival = pd.DataFrame()
for i in range(len(features_survival)):
    data_survival[features_survival[i]] = data[features_survival[i]]
data_cellcycle = pd.DataFrame()
for i in range(len(features_cellcycle)):
    data_cellcycle[features_cellcycle[i]] = data[features_cellcycle[i]]

data_survival_flt = data_survival[data_survival['mean_n_filtered'] >= 200].copy().reset_index(drop=True)
data_survival_flt = data_survival_flt.sort_values(by='log2fc_pos_vs_neg', ascending=False).reset_index(drop=True)
data_survival_flt.index = data_survival_flt['label']
data_survival_flt = data_survival_flt.drop(columns=['mean_n_filtered', 'label'])

# heat map
fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_survival_flt, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')
plt.savefig('%s/%s_pos_vs_neg_heatmap.pdf' % (output_dir, batch))
plt.show()