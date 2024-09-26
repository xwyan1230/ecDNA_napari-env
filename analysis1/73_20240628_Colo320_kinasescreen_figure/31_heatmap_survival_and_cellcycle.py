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
"""exclude = ['Barasertib', 'Alisertib', 'AZD2811', 'PF-03814735', 'GSK269962A', 'Y-33075', 'BDP5290', 'Midostaurin', 'FF-10101',
           'MRT67307', 'Axitinib', 'XMU-MP-1', 'XL413']"""
exclude = []
exclude_index = [78, 85, 134, 154, 188]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
label_seq = pd.read_csv('%s/cluster_label_hoechst_survival_only_exclude.txt' % (output_dir), na_values=['.'], sep='\t')

df = pd.DataFrame()
for batch in batches:
    features = ['mean_hoechst']
    data = pd.read_csv('%s/%s/%s_average.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    df['treatment'] = data['treatment']
    for feature in features:
        df['%s_%s' % (batch, feature)] = data[feature]

for batch in batches:
    features = ['mean_per_n_G1', 'mean_per_n_S', 'mean_per_n_G2M']
    data = pd.read_csv('%s/%s/%s_cellcycle_average.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    df['treatment'] = data['treatment']
    for feature in features:
        df['%s_%s' % (batch, feature)] = data[feature]

df_flt = df[~df['treatment'].isin(exclude)].copy().reset_index(drop=True)
print(len(df_flt))
df = df_flt.drop(exclude_index).reset_index(drop=True)
df_flt = df.copy()
print(len(df))

df_flt.index = df_flt['treatment']
df_flt = df_flt.drop('treatment', axis=1)

scaler = StandardScaler()
# data_heatmap_scale = pd.DataFrame(scaler.transform(df_flt))
data_heatmap_scale = pd.DataFrame(scaler.fit_transform(df_flt))
print(data_heatmap_scale.shape)
data_heatmap_scale.columns = df_flt.columns

# hierarchical Clustering
# DHC
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.

plt.figure(figsize=(20, 3))
clusters = shc.linkage(data_heatmap_scale, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/ahc_survival_and_cellcycle.pdf' % output_dir)
plt.show()
nodes = R['ivl']

df_sort = pd.DataFrame(columns=df_flt.columns)
label_sort = []
for i in range(len(df_flt)):
    df_sort.loc[len(df_sort.index)] = \
        df_flt.iloc[int(dat.list_invert(nodes)[i])]
    label_sort.append(df['treatment'].tolist()[int(dat.list_invert(nodes)[i])])

df_sort.index = label_sort

df_heatmap_sort = pd.DataFrame(columns=data_heatmap_scale.columns)
label_sort = []
for i in range(len(data_heatmap_scale)):
    df_heatmap_sort.loc[len(df_heatmap_sort.index)] = \
        data_heatmap_scale.iloc[int(dat.list_invert(nodes)[i])]
    label_sort.append(df['treatment'].tolist()[int(dat.list_invert(nodes)[i])])

df_heatmap_sort.index = label_sort

pd_label = pd.DataFrame()
pd_label['node'] = dat.list_invert(nodes)
pd_label['label'] = label_sort
pd_label['target'] = ['DMSO' if pd_label['label'][i] == 'DMSO'else target_df[target_df['compound name'] == pd_label['label'][i]]['TargetGene'].tolist()[0] for i in range(len(pd_label))]
pd_label_heatmap = pd_label.copy()

pd_label.to_csv('%s/cluster_label_survival_and_cellcycle.txt' % (output_dir), index=False, sep='\t')

# heat map
fig, ax = plt.subplots(figsize=(6, 15))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_heatmap_scale, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')
# plt.savefig('%s/heatmap.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_heatmap_sort, cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm', yticklabels=False)  #  yticklabels=False
plt.savefig('%s/heatmap_sort_survival_and_cellcycle.pdf' % output_dir)
plt.show()

"""fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_heatmap_sort[0:34], cbar=0, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_sort_enlarge.pdf' % output_dir)
plt.show()"""

"""pd_ctrl = pd.DataFrame()
pd_ctrl['label'] = label_sort
pd_ctrl['ctrl'] = [1 if pd_ctrl['label'][i] == 'DMSO'else 0 for i in range(len(pd_ctrl))]
pd_ctrl.index = pd_ctrl['label']
pd_ctrl = pd_ctrl.drop('label', axis=1)

norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), "darkblue"],
          [norm(-0.6), "lightgrey"],
          [norm( 0.6), "lightgrey"],
          [norm( 1.0), rainboo_colors[6]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(pd_ctrl, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap, yticklabels=False)  #  yticklabels=False
plt.savefig('%s/heatmap_DMSO.pdf' % output_dir)
plt.show()

gene_lst = ['PLK1', 'MTOR', 'PTK2', 'CDK2', 'CDK4', 'CDK7', 'CDK9', 'CDK12']
pd_genes = pd.DataFrame()
pd_genes['target'] = pd_label['target']
for gene in gene_lst:
    pd_genes[gene] = [1 if gene in str(pd_genes['target'][i]) else 0 for i in range(len(pd_genes))]
pd_genes.index = pd_label['label']
pd_genes = pd_genes.drop('target', axis=1)

fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(pd_genes[0:34], cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_sort_enlarge_gene.pdf' % output_dir)
plt.show()"""