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

batchA = 'point5uM_24hr'
batchB = 'point5uM_48hr'
batchC = '5uM_24hr'
batchD = '5uM_48hr'
features_heatmap = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'log2fc_per_G1', 'log2fc_per_S', 'log2fc_per_G2M']
features_plot = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis']
features = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S', 'log2fc_per_G2M', 'log2fc_cytokinesis']
features_plot_survival = ['log2fc_fov_hoechst', 'log2fc_n_filtered']
features_plot_survival = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'log2fc_n_neg', 'log2fc_n_pos', 'log2fc_pos_vs_neg']
features_plot_cellcycle = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis', 'cc_score_neg', 'log2fc_per_neg_G1', 'log2fc_per_neg_G1S',
                           'log2fc_per_neg_S', 'log2fc_per_neg_G2M', 'cc_score_pos',
                           'log2fc_per_pos_G1', 'log2fc_per_pos_G1S', 'log2fc_per_pos_S', 'log2fc_per_pos_G2M']

features_plot = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis']


def data_load(data_dir, batch):
    data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
    data = pd.concat([data, data_cc], axis=1)
    data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
    data_heatmap = pd.DataFrame()
    for i in range(len(features_heatmap)):
        data_heatmap['%s_%s' % (features_heatmap[i], batch)] = data[features_heatmap[i]]
    data_heatmap_flt = data_heatmap.copy()
    for feature in features:
        if '%s_%s' % (feature, batch) in data_heatmap.columns:
            data_heatmap_flt.loc[(data['rep'] == 2) & (data['mean_n_filtered'] < 200), '%s_%s' % (feature, batch)] = np.nan
    data_plot = pd.DataFrame()
    for i in range(len(features_plot)):
        data_plot['%s_%s' % (features_plot[i], batch)] = data[features_plot[i]]
    data_plot_flt = data_plot.copy()
    for feature in features:
        if '%s_%s' % (feature, batch) in data_plot.columns:
            data_plot_flt.loc[(data['rep'] == 2) & (data['mean_n_filtered'] < 200), '%s_%s' % (feature, batch)] = np.nan
    data_plot_survival = pd.DataFrame()
    for i in range(len(features_plot_survival)):
        data_plot_survival['%s_%s' % (features_plot_survival[i], batch)] = data[features_plot_survival[i]]
    data_plot_cellcycle = pd.DataFrame()
    for i in range(len(features_plot_cellcycle)):
        data_plot_cellcycle['%s_%s' % (features_plot_cellcycle[i], batch)] = data[features_plot_cellcycle[i]]
    return data, data_heatmap, data_heatmap_flt, data_plot, data_plot_flt, data_plot_survival, data_plot_cellcycle


dataA, dataA_heatmap, dataA_heatmap_flt, dataA_plot, dataA_plot_flt, dataA_plot_survival, dataA_plot_cellcycle = data_load(data_dir, batchA)
dataB, dataB_heatmap, dataB_heatmap_flt, dataB_plot, dataB_plot_flt, dataB_plot_survival, dataB_plot_cellcycle = data_load(data_dir, batchB)
dataC, dataC_heatmap, dataC_heatmap_flt, dataC_plot, dataC_plot_flt, dataC_plot_survival, dataC_plot_cellcycle = data_load(data_dir, batchC)
dataD, dataD_heatmap, dataD_heatmap_flt, dataD_plot, dataD_plot_flt, dataD_plot_survival, dataD_plot_cellcycle = data_load(data_dir, batchD)

data_heatmap = pd.concat([dataA_heatmap, dataB_heatmap, dataC_heatmap, dataD_heatmap], axis=1)
data_plot = pd.concat([dataA_plot, dataB_plot, dataC_plot, dataD_plot], axis=1)
data_plot_flt = pd.concat([dataA_plot_flt, dataB_plot_flt, dataC_plot_flt, dataD_plot_flt], axis=1)
data_plot_survival = pd.concat([dataA_plot_survival, dataB_plot_survival, dataC_plot_survival, dataD_plot_survival], axis=1)
data_plot_cellcycle = pd.concat([dataA_plot_cellcycle, dataB_plot_cellcycle, dataC_plot_cellcycle, dataD_plot_cellcycle], axis=1)
print(len(data_heatmap))

cmap = pcm.get_cmap('Spectral')
line_colors = []
for i in np.arange(0, 1, 1/len(dataA)):
    line_colors.append(cmap(i))
line_colors.append((0.30, 0.30, 0.30, 1.0))

"""scaler = StandardScaler()
data_heatmap_scale = pd.DataFrame(scaler.transform(data_heatmap))
data_heatmap_scale.columns = data_heatmap.columns"""

# hierarchical Clustering
# DHC
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.

plt.figure(figsize=(20, 3))
clusters = shc.linkage(data_heatmap, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/ahc.pdf' % output_dir)
plt.show()
nodes = R['ivl']

data_plot_sort = pd.DataFrame(columns=data_plot.columns)
data_plot_survival_sort = pd.DataFrame(columns=data_plot_survival.columns)
data_plot_cellcycle_sort = pd.DataFrame(columns=data_plot_cellcycle.columns)
label_sort = []
target_sort = []
for i in range(len(data_plot)):
    data_plot_sort.loc[len(data_plot_sort.index)] = \
        data_plot.iloc[int(dat.list_invert(nodes)[i])]
    data_plot_survival_sort.loc[len(data_plot_survival_sort.index)] = \
        data_plot_survival.iloc[int(dat.list_invert(nodes)[i])]
    data_plot_cellcycle_sort.loc[len(data_plot_cellcycle_sort.index)] = \
        data_plot_cellcycle.iloc[int(dat.list_invert(nodes)[i])]
    label_sort.append(dataA['label'].tolist()[int(dat.list_invert(nodes)[i])])
    target_sort.append(dataA['target'].tolist()[int(dat.list_invert(nodes)[i])])

data_plot.index = dataA['label'].tolist()

pd_label = pd.DataFrame()
pd_label['node'] = dat.list_invert(nodes)
pd_label['label'] = label_sort
pd_label['target'] = target_sort
pd_label_heatmap = pd_label.copy()
"""search_lst = ['DMSO', 'PLK1', 'MELK', 'CDK9', 'CDK12|CDK13', 'AKT|Akt', 'HPK1', 'CDK7',
              'CDK4|CDK6', 'NEK2|Nek2', 'CDK2', 'mTOR', 'CHK1', 'BTK', 'ALK', 'FLT3', 'FAK',
              'HER2|ERBB2|ErbB-2|ErbB2', 'Axl', 'Aurora', 'ROCK']"""
search_lst = ['ATM', 'ATR']
for search in search_lst:
    pd_label_heatmap[search] = [0] * len(pd_label_heatmap)
    pd_label_heatmap.loc[pd_label_heatmap['target'].str.contains(search), search] = 1
pd_label_heatmap = pd_label_heatmap.drop(['node', 'label', 'target'], axis=1)
print(pd_label_heatmap.head())

# pd_label.to_csv('%s/cluster_label.txt' % (output_dir), index=False, sep='\t')

# heat map
fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_plot, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')
# plt.savefig('%s/heatmap.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_plot_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm', yticklabels=False)
# plt.savefig('%s/heatmap_sort.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_plot_survival_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm',
                  yticklabels=False)
# plt.savefig('%s/heatmap_survival_sort_all.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_plot_cellcycle_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm',
                  yticklabels=False)
# plt.savefig('%s/heatmap_cellcycle_sort_all.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(pd_label_heatmap, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap='coolwarm',
                  yticklabels=False)
plt.savefig('%s/heatmap_label_ATM_ATR.pdf' % output_dir)
plt.show()



print("DONE!")