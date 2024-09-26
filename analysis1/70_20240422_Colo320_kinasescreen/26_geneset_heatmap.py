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

batchA = 'point5uM_24hr'
batchB = 'point5uM_48hr'
batchC = '5uM_24hr'
batchD = '5uM_48hr'

features_heatmap = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis']
# features_survival = ['log2fc_fov_hoechst', 'log2fc_n_filtered']
# features_cellcycle = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S', 'log2fc_per_G2M', 'log2fc_cytokinesis']

features_survival = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'log2fc_n_neg', 'log2fc_n_pos', 'log2fc_pos_vs_neg']
features_cellcycle = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
                'log2fc_per_G2M', 'log2fc_cytokinesis', 'cc_score_neg', 'log2fc_per_neg_G1', 'log2fc_per_neg_G1S',
                           'log2fc_per_neg_S', 'log2fc_per_neg_G2M', 'cc_score_pos',
                           'log2fc_per_pos_G1', 'log2fc_per_pos_G1S', 'log2fc_per_pos_S', 'log2fc_per_pos_G2M']
# point5uM 24hr survival
# search_lst = ['PLK1', 'CDK13', 'MELK', 'CDK9', 'BRD4', 'AKT1', 'CDK12', 'PIK3CA']
# point5uM 48hr survival
# search_lst = ['BRD4', 'CDK12', 'AKT1', 'CDK9', 'MELK', 'CDK13', 'PLK1', 'PIK3CA']
# 5uM 24hr survival
# search_lst = ['BRD4', 'SRC', 'FLT3', 'PLK1', 'PTK2', 'CDK12', 'CDK9', 'CDK2', 'PDGFRB', 'ALK', 'AURKB', 'FGFR1', 'MTOR', 'MAP4K1', 'IGF1R', 'MAP2K7', 'CDK4', 'SIK1', 'CDK6', 'SIK3', 'MELK', 'NEK2', 'SIK2', 'ABL1', 'PIK3CA', 'ABL2', 'CDK13', 'BTK', 'PDGFRA', 'AKT1', 'CDK7', 'AURKA', 'KDR', 'INSR']
# 5uM 48hr survival
# search_lst = ['MTOR', 'AURKA', 'IGF1R', 'PDGFRA', 'ABL2', 'ALK', 'CDK9', 'AXL', 'PIK3CA', 'EGFR', 'INSR', 'CDK2', 'CHK1', 'MAP2K7', 'FLT3', 'ERBB2', 'SIK1', 'AKT1', 'BTK', 'AURKB', 'PTK2', 'NEK2', 'CDK6', 'SIK3', 'FGFR1', 'CDK12', 'KDR', 'BRD4', 'CDK7', 'MELK', 'CDK4', 'SIK2', 'ABL1', 'PLK1', 'CDK13', 'MAP4K1', 'PDGFRB', 'SRC']
# survival_hit_14
# search_lst = ['PLK1', 'MTOR', 'CHK1', 'CDK9', 'CDK12', 'CDK13', 'CDK4', 'CDK6', 'CDK7', 'AURKA', 'AURKB', 'IGF1R', 'INSR', 'PTK2']
# 5uM 24hr pos_vs_neg
# search_lst = ['ALK', 'CDK4', 'IGF1R', 'PDGFRB', 'ABL1', 'MAP3K5', 'ABL2', 'KDR', 'ROCK2', 'FGFR1', 'PTK2', 'CDK6', 'SRC', 'MTOR', 'PDGFRA', 'FLT3', 'INSR', 'ROCK1']
# 5uM 48hr pos_vs_neg
# search_lst = ['AURKB', 'PTK2', 'FLT4', 'ABL2', 'ROCK1', 'SRC', 'KDR', 'CDK13', 'HCK', 'IGF1R', 'LRRK2', 'CDK9', 'CDK12', 'AURKA', 'FLT3', 'PDGFRA', 'KIT', 'FGFR1', 'MTOR', 'BCR', 'ABL1', 'CDK7', 'BRD4', 'CSF1R', 'ALK', 'MAP3K5', 'PDGFRB', 'CDK2', 'ROCK2', 'PLK1', 'INSR', 'MAP4K1', 'FLT1', 'LYN', 'CHK1', 'FGR']
# 5uM 48hr pos_vs_neg_reverse
# search_lst = ['AKT1', 'JAK1', 'CDK13', 'JAK2', 'CDK2', 'JAK3', 'CDK12', 'GSK3A', 'GSK3B']
# 5uM 24hr pos_vs_neg_reverse
# search_lst = ['CDK2', 'PIK3CA', 'JAK3', 'CDK12', 'CDK9', 'CDK13', 'JAK1', 'CDK7', 'JAK2']
# DDR
search_lst = ['ATM', 'ATR', 'CHK1', 'CHK2', 'WEE1']
search_lst.sort()
print(search_lst)
hit_name = 'DDR'


def data_load(data_dir, batch):
    data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
    data = pd.concat([data, data_cc], axis=1)
    data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
    count_lst = []
    data_flt = pd.DataFrame()
    for search in search_lst:
        data_temp = data[data['target'].str.contains(search)].copy().reset_index(drop=True)
        count_lst.append(len(data_temp))
        data_temp.loc[len(data_temp.index)] = [np.nan] * len(data_temp.columns)
        data_flt = pd.concat([data_flt, data_temp], axis=0)
    data_heatmap = pd.DataFrame()
    for i in range(len(features_heatmap)):
        data_heatmap['%s_%s' % (features_heatmap[i], batch)] = data_flt[features_heatmap[i]]
    data_survival = pd.DataFrame()
    for i in range(len(features_survival)):
        data_survival['%s_%s' % (features_survival[i], batch)] = data_flt[features_survival[i]]
    data_cellcycle = pd.DataFrame()
    for i in range(len(features_cellcycle)):
        data_cellcycle['%s_%s' % (features_cellcycle[i], batch)] = data_flt[features_cellcycle[i]]
    return data, data_heatmap, count_lst, data_survival, data_cellcycle


dataA, dataA_heatmap, countA, dataA_survival, dataA_cellcycle = data_load(data_dir, batchA)
dataB, dataB_heatmap, _, dataB_survival, dataB_cellcycle = data_load(data_dir, batchB)
dataC, dataC_heatmap, _, dataC_survival, dataC_cellcycle = data_load(data_dir, batchC)
dataD, dataD_heatmap, _, dataD_survival, dataD_cellcycle = data_load(data_dir, batchD)

data_heatmap = pd.concat([dataA_heatmap, dataB_heatmap, dataC_heatmap, dataD_heatmap], axis=1)
data_survival = pd.concat([dataA_survival, dataB_survival, dataC_survival, dataD_survival], axis=1)
data_cellcycle = pd.concat([dataA_cellcycle, dataB_cellcycle, dataC_cellcycle, dataD_cellcycle], axis=1)
y_axis_labels = []
for i in range(len(search_lst)):
    y_axis_labels = y_axis_labels + [search_lst[i]] + [''] * countA[i]

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm',
                  yticklabels=y_axis_labels)
plt.savefig('%s/%s_heatmap.pdf' % (output_dir, hit_name))
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_survival, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm',
                  yticklabels=y_axis_labels)
plt.savefig('%s/%s_heatmap_survival.pdf' % (output_dir, hit_name))
plt.show()

fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_cellcycle, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm',
                  yticklabels=y_axis_labels)
plt.savefig('%s/%s_heatmap_cellcycle.pdf' % (output_dir, hit_name))
plt.show()


