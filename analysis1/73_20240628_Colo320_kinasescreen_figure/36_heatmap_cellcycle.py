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
exclude_index = [154, 188]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), na_values=['.'], sep='\t')

df = pd.DataFrame()
features = ['delta_n_G1', 'delta_n_G1S', 'delta_n_S', 'delta_n_G2M', 'delta_n_G2MG1']
for batch in batches:
    data = pd.read_csv('%s/%s/%s_cellcycle_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
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

"""# hierarchical Clustering
# DHC
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.
plt.figure(figsize=(20, 3))
clusters = shc.linkage(df_num_cluster, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/ahc_rescale_survival_hoechst.pdf' % output_dir)
plt.show()
nodes = R['ivl']

df_label_sort = pd.DataFrame(columns=df_label.columns)
for i in range(len(df_label)):
    df_label_sort.loc[len(df_label_sort.index)] = df_label.iloc[int(dat.list_invert(nodes)[i])]
df_label_sort['seq1'] = df_label_sort.index
seq = []
for i in range(len(df_label_sort)):
    if df_label_sort['seq1'][i] >= (n_len-17):
        seq.append(i-(n_len-17))
    elif df_label_sort['seq1'][i] <= 35:
        seq.append(52-i)
    else:
        seq.append(i+17)
df_label_sort['seq'] = seq
df_label_sort = df_label_sort.sort_values(by='index').reset_index(drop=True)
df_label_sort['target'] = ['DMSO' if df_label_sort['treatment'][i] == 'DMSO'else target_df[target_df['compound name'] == df_label_sort['treatment'][i]]['TargetGene'].tolist()[0] for i in range(len(df_label_sort))]
print(df_label_sort)

df_cluster = df_num_cluster.copy()
df_cluster['index'] = df_label_sort['index']
df_cluster['treatment'] = df_label_sort['treatment']
df_cluster['seq'] = df_label_sort['seq']
df_cluster['target'] = df_label_sort['target']
df_cluster.to_csv('%s/average_survival_hoechst_scale.txt' % (output_dir), index=False, sep='\t')

df_cluster_sort = df_cluster.copy()
df_cluster_sort = df_cluster_sort.sort_values(by='seq').reset_index(drop=True)
df_cluster_sort.to_csv('%s/average_survival_hoechst_scale_sort.txt' % (output_dir), index=False, sep='\t')

df_num_cluster_sort = df_num_cluster.copy()
df_num_cluster_sort['treatment'] = df_label_sort['treatment'].tolist()
df_num_cluster_sort['seq'] = df_label_sort['seq'].tolist()
df_num_cluster_sort = df_num_cluster_sort.sort_values(by='seq')
df_num_cluster_sort.index = df_num_cluster_sort['treatment']
df_num_cluster_sort = df_num_cluster_sort.drop(['seq', 'treatment'], axis=1)"""

# heat map
fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_num_cluster_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_cellcycle_survival_seq.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_num_cluster_sort[0:50], cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_cellcycle_survival_seq_enlarge.pdf' % output_dir)
plt.show()