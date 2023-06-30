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

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder
if not os.path.exists("%s/cluster/" % output_dir):
    os.makedirs("%s/cluster/" % output_dir)

feature_lst = ['mean_area_nuclear', 'mean_averageD_point2', 'std_averageD_point2',
               'mean_mean_int_DNAFISH', 'mean_n_ecDNA', 'mean_r10', 'mean_r16']

exclude = ['D9']
features_heatmap = feature_lst
df_mean = pd.read_csv('%ssummary_mean.txt' % data_dir, na_values=['.'], sep='\t')
data_mean = df_mean[~df_mean['sample'].isin(exclude)].copy().reset_index(drop=True)

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

plt.figure(figsize=(10, 3))
clusters = shc.linkage(data_mean_feature_scale, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/cluster/ahc.pdf' % output_dir)
plt.close()
nodes = R['ivl']

data_mean_feature_scale_sort = pd.DataFrame(columns=data_mean_feature_scale.columns)
for i in range(len(data_mean_feature_scale)):
    data_mean_feature_scale_sort.loc[len(data_mean_feature_scale_sort.index)] = \
        data_mean_feature_scale.iloc[int(dat.list_invert(nodes)[i])]

# heat map
plt.subplots(figsize=(8, 15))
ax1 = sns.heatmap(data_mean_feature_scale, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis')
plt.savefig('%s/cluster/heatmap.pdf' % output_dir)
plt.close()

plt.subplots(figsize=(8, 15))
ax1 = sns.heatmap(data_mean_feature_scale_sort, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis',
                  yticklabels=False)
plt.savefig('%s/cluster/heatmap_sort.pdf' % output_dir)
plt.close()

data_feature = data_mean.copy()
all_feature = data_mean.columns
drop_feature = [i for i in all_feature if i not in feature_lst]
data_feature = data_feature.drop(drop_feature, axis=1)

scaler = StandardScaler()
data_feature_scale = pd.DataFrame(scaler.fit_transform(data_feature))
data_feature_scale.columns = data_feature.columns

# pca
# 1. rescale the data to a [0,1] range
# 2. standardize the data to have a zero mean and unit standard deviation (currently using this one)
# data_feature_scale = (data_feature - data_feature.min())/(data_feature.max() - data_feature.min())
plt.subplots(figsize=(12, 10))
pca = sklearnPCA(n_components=2)
transformed = pd.DataFrame(pca.fit_transform(data_feature_scale))
eigenvalues = pca.explained_variance_
print(eigenvalues)

for i in range(len(data_mean)):
    genelist = data_mean['gene'].tolist()
    if genelist[i] != 'Ctrl':
        plt.scatter(transformed[data_mean['gene'] == genelist[i]][0], transformed[data_mean['gene'] == genelist[i]][1],
                    label=genelist[i], color=line_colors[i], alpha=0.5, s=5)
    else:
        plt.scatter(transformed[data_mean['gene'] == genelist[i]][0], transformed[data_mean['gene'] == genelist[i]][1],
                    label=genelist[i], color=line_colors[-1], alpha=0.5, s=5)
plt.legend()
plt.savefig('%s/cluster/pca.pdf' % output_dir)
plt.show()

# lda
plt.subplots(figsize=(12, 10))
lda = LDA(n_components=2)
lda_transformed = pd.DataFrame(lda.fit_transform(data_feature_scale, data_mean['gene']))

# Plot all three series
for i in range(len(genelist)):
    if genelist[i] != 'Ctrl':
        plt.scatter(lda_transformed[data_mean['gene'] == genelist[i]][0], lda_transformed[data_mean['gene'] == genelist[i]][1],
                    label=genelist[i], c=line_colors[i], alpha=0.5, s=5)
    else:
        plt.scatter(lda_transformed[data_mean['gene'] == genelist[i]][0], lda_transformed[data_mean['gene'] == genelist[i]][1],
                    label=genelist[i], c=line_colors[-1], alpha=0.5, s=5)

plt.legend()
plt.savefig('%s/cluster/lda.pdf' % output_dir)
plt.show()

"""# parallel coordinates
print("Plotting pc...")
# Select features to include in the plot
plot_feat = data_feature.columns
# Concat classes with the normalized data
data_norm = pd.concat([data_feature_scale[plot_feat], data['sample']], axis=1)
# Perform parallel coordinate plot
plt.subplots(figsize=(20, 4))
parallel_coordinates(data_norm, 'sample', color=line_color_4, alpha=0.3)
# Display legend and show plot
plt.legend(loc=3)
plt.savefig('%s/pc.pdf' % save_folder)
plt.close()"""

print("DONE!")