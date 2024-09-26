import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import matplotlib
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
exclude_index = [78, 85, 134, 154, 188]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale.txt' % (output_dir), na_values=['.'], sep='\t')

df = pd.DataFrame()
features = ['log2_mean_hoechst']
for batch in batches:
    data = pd.read_csv('%s/%s/%s_average.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data_temp = pd.DataFrame()
    data_temp['index'] = data.index
    data_temp['treatment'] = data['treatment']
    for feature in features:
        data_temp[feature] = data[feature]
    df = pd.concat([df, data_temp], axis=0).reset_index(drop=True)

df_label = data_temp.copy()
df_label = df_label.drop(features, axis=1)
df_label['survival_seq'] = survival_seq['seq']
df_label_drop = df_label.copy()
df_label_drop = df_label_drop.drop(exclude_index).reset_index(drop=True)
print(len(df_label_drop))

df_num = df.copy()
df_num = df_num[~df_num['index'].isin(exclude_index)].copy().reset_index(drop=True)
df_num = df_num.drop(['index', 'treatment'], axis=1)
scaler = StandardScaler()
df_num_scale = pd.DataFrame(scaler.fit_transform(df_num))
df_num_scale.columns = df_num.columns
df['%s_scale' % feature] = df_num_scale[feature]

n_len = int(len(df_num)/4)
print(n_len)
df_num_cluster = pd.DataFrame()
for feature in features:
    for i in range(4):
        df_num_cluster['%s_%s' % (feature, i)] = df_num_scale.iloc[i*n_len:(i+1)*n_len][feature].reset_index(drop=True)
df_num_cluster['survival_seq'] = df_label_drop['survival_seq']
df_num_cluster['treatment'] = df_label_drop['treatment']
print(df_num_cluster)
print(len(df_num_cluster))

df_num_cluster_sort = df_num_cluster.copy()
df_num_cluster_sort = df_num_cluster_sort.sort_values(by='survival_seq')
df_num_cluster_sort.index = df_num_cluster_sort['treatment']
print(df_num_cluster_sort)
df_num_cluster_sort = df_num_cluster_sort.drop(['survival_seq', 'treatment'], axis=1)

# heat map
fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_num_cluster_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_hoechst_survival_seq.pdf' % output_dir)
plt.show()

pd_ctrl = df_label_drop.copy()
pd_ctrl['ctrl'] = [1 if pd_ctrl['treatment'][i] == 'DMSO'else 0 for i in range(len(pd_ctrl))]
pd_ctrl = pd_ctrl.sort_values(by='survival_seq')
pd_ctrl = pd_ctrl.drop(['index', 'treatment', 'survival_seq'], axis=1)

norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), "darkblue"],
          [norm(-0.6), "#d0d2d3"],
          [norm( 0.6), "#d0d2d3"],
          [norm( 1.0), rainboo_colors[6]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(pd_ctrl, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap, yticklabels=False)  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_hoechst_survival_seq_DMSO.pdf' % output_dir)
plt.show()