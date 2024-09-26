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


def data_load(data_dir, batch):
    data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
    data = pd.concat([data, data_cc], axis=1)
    data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
    data['cc_ratio'] = (data['cc_score_pos'] + 1) / (data['cc_score_neg'] + 1)
    data_heatmap = pd.DataFrame()
    data_heatmap['log2fc_fov_hoechst_%s' % batch] = [1 if data['log2fc_fov_hoechst'][i] < -1 else 0 for i in range(len(data))]
    # data_heatmap['log2fc_n_filtered_%s' % batch] = [1 if data['log2fc_n_filtered'][i] < -1 else 0 for i in range(len(data))]
    # data_heatmap['cc_score_%s' % batch] = [1 if data['cc_score'][i] > 1 else 0 for i in range(len(data))]
    data_heatmap['log2fc_pos_vs_neg_%s' % batch] = [1 if (data['mean_n_filtered'][i] > 200) & ((data['log2fc_pos_vs_neg'][i] < -0.6)|(data['log2fc_pos_vs_neg'][i] > 0.6)) else 0 for i in
                                                     range(len(data))]
    # data_heatmap['cc_ratio_%s' % batch] = [1 if (data['mean_n_filtered'][i] > 200) & ((data['cc_ratio'][i] < 0.5) | (data['cc_ratio'][i] > 1.5)) else 0 for i in range(len(data))]
    return data, data_heatmap


dataA, dataA_heatmap = data_load(data_dir, batchA)
dataB, dataB_heatmap = data_load(data_dir, batchB)
dataC, dataC_heatmap = data_load(data_dir, batchC)
dataD, dataD_heatmap = data_load(data_dir, batchD)

data_heatmap = pd.concat([dataA_heatmap, dataB_heatmap, dataC_heatmap, dataD_heatmap], axis=1)
data_heatmap['sum'] = data_heatmap.sum(axis=1)
data_heatmap.index = dataA['label']
data_heatmap = data_heatmap.sort_values(by='sum')
print(data_heatmap.head())

cmap = pcm.get_cmap('Spectral')
line_colors = []
for i in np.arange(0, 1, 1/len(dataA)):
    line_colors.append(cmap(i))
line_colors.append((0.30, 0.30, 0.30, 1.0))

# heat map
fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_heatmap, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap='coolwarm')
plt.savefig('%s/heatmap_hit.pdf' % output_dir)
plt.show()

print(list(reversed(data_heatmap.index)))

print("DONE!")