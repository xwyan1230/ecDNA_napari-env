import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import skimage.io as skio
import shared.display as dis
import shared.dataframe as dat
import shared.image as ima
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import napari
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220826_neuroblastoma/"
group = 'ABC'
group_lst = ['A', 'B', 'C']
save_path = '%sv1_img/%s/' % (master_folder, group)
if not os.path.exists(save_path):
    os.makedirs(save_path)
version = 1
cmap = matplotlib.cm.get_cmap('Spectral')
features_heatmap = ['area_nuclear', 'n_ecDNA', 'total_area_ratio_ecDNA', 'relative_r_area', 'dis_to_hub_area',
                    'percentage_area_ratio_n_half', 'cum_area_ratio_n_half']

data_mean = pd.read_csv('%s%s_v%s_mean.txt' % (master_folder, group, version), na_values=['.'], sep='\t')

cmap = pcm.get_cmap('Spectral')
line_colors = []
for i in np.arange(0, 1, 1/len(data_mean)):
    line_colors.append(cmap(i))
line_colors.append((0.30, 0.30, 0.30, 1.0))

data_mean_feature = data_mean.copy()
all_mean_feature = data_mean.columns
drop_feature = [i for i in all_mean_feature if i not in features_heatmap]
data_mean_feature = data_mean_feature.drop(drop_feature, axis=1)

scaler = StandardScaler()
data_mean_feature_scale = pd.DataFrame(scaler.fit_transform(data_mean_feature))
data_mean_feature_scale.columns = data_mean_feature.columns

# hierarchical Clustering
# AHC (agglomerative, bottom-up)
# the other way is DHC (divisive, top-down). DHC works better when you have fewer but larger clusters, hence it's more
# computationally expensive. AHC is fitted for when you have many smaller clusters. It is computationally simpler, more
# used and more available.

plt.figure(figsize=(60, 8))
clusters = shc.linkage(data_mean_feature_scale, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/ahc.pdf' % save_path)
plt.close()
nodes = R['ivl']
sample_nodes = [data_mean['sample'].tolist()[int(i)] for i in nodes]
group_nodes = [data_mean['group'].tolist()[int(i)] for i in nodes]
group_nodes_number = []
for i in group_nodes:
    if i == 'A':
        group_nodes_number.append(-5)
    elif i == 'B':
        group_nodes_number.append(0)
    elif i == 'C':
        group_nodes_number.append(5)
print(sample_nodes)
print(group_nodes)

data_mean_feature_scale_sort = pd.DataFrame(columns=data_mean_feature_scale.columns)
for i in range(len(data_mean_feature_scale)):
    data_mean_feature_scale_sort.loc[len(data_mean_feature_scale_sort.index)] = \
        data_mean_feature_scale.iloc[int(dat.list_invert(nodes)[i])]
data_mean_feature_scale_sort['group'] = group_nodes_number

# heat map
plt.subplots(figsize=(8, 60))
ax1 = sns.heatmap(data_mean_feature_scale, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis')
plt.savefig('%s/heatmap.pdf' % save_path)
plt.close()

plt.subplots(figsize=(8, 60))
ax1 = sns.heatmap(data_mean_feature_scale_sort, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis',
                  yticklabels=False)
plt.savefig('%s/heatmap_sort.pdf' % save_path)
plt.close()

data_feature = data_mean.copy()
all_feature = data_mean.columns
drop_feature = [i for i in all_feature if i not in features_heatmap]
data_feature = data_feature.drop(drop_feature, axis=1)

scaler = StandardScaler()
data_feature_scale = pd.DataFrame(scaler.fit_transform(data_feature))
data_feature_scale.columns = data_feature.columns

line_colors = []
for i in np.arange(0, 1, 1/len(group_lst)):
    line_colors.append(cmap(i))
line_colors.append((0.30, 0.30, 0.30, 1.0))

# pca
# 1. rescale the data to a [0,1] range
# 2. standardize the data to have a zero mean and unit standard deviation (currently using this one)
# data_feature_scale = (data_feature - data_feature.min())/(data_feature.max() - data_feature.min())
plt.subplots(figsize=(12, 10))
pca = sklearnPCA(n_components=2)
transformed = pd.DataFrame(pca.fit_transform(data_feature_scale))
eigenvalues = pca.explained_variance_
print(eigenvalues)

for i in range(len(group_lst)):
    plt.scatter(transformed[data_mean['group'] == group_lst[i]][0], transformed[data_mean['group'] == group_lst[i]][1],
                label=group_lst[i], color=line_colors[i], alpha=0.8, s=7)

plt.legend()
plt.savefig('%s/pca.pdf' % save_path)
plt.close()

# lda
plt.subplots(figsize=(12, 10))
lda = LDA(n_components=2)
lda_transformed = pd.DataFrame(lda.fit_transform(data_feature_scale, data_mean['group']))

# Plot all three series
for i in range(len(group_lst)):
    plt.scatter(lda_transformed[data_mean['group'] == group_lst[i]][0], lda_transformed[data_mean['group'] == group_lst[i]][1],
                label=group_lst[i], c=line_colors[i], alpha=0.8, s=7)

plt.legend()
plt.savefig('%s/lda.pdf' % save_path)
plt.close()

print("DONE!")
