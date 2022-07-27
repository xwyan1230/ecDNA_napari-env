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

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220526_flowFISH_topHits_screen/"
features_heatmap = ['area_nuclear', 'mean_int_nuclear', 'total_int_nuclear', 'mean_int_IF', 'total_int_IF', 'n_ecDNA',
                    'mean_int_DNAFISH_norm', 'total_int_DNAFISH_norm', 'mean_int_ecDNA_norm', 'total_int_ecDNA_norm',
                    'mean_int_ind_ecDNA_norm', 'total_int_ind_ecDNA_norm', 'area_ratio_ind_ecDNA',
                    'area_ratio_ecDNA', 'max_area_ratio_ecDNA', 'radial_center', 'cum_area_ratio_n_half',
                    'cum_int_norm_n_half', 'dis_to_hub_area', 'dis_to_hub_int_norm', 'g_value', 'angle_value', ]
features = ['area_nuclear', 'mean_int_DNAFISH_norm', 'mean_int_IF', 'radial_center', 'area_ratio_ecDNA',
            'mean_int_ecDNA_norm', 'max_area_ratio_ecDNA', 'n_ecDNA', 'cum_area_ratio_n_half',
            'cum_int_norm_n_half', 'dis_to_hub_int_norm', 'g_value', 'angle_value']
save_folder = "%scluster/" % master_folder
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

data = pd.read_csv(("%ssummary/summary.txt" % master_folder), na_values=['.'], sep='\t')
data_mean = pd.read_csv(("%ssummary/summary_mean.txt" % master_folder), na_values=['.'], sep='\t')
data_gene = pd.read_csv("%sgene.txt" % master_folder, na_values=['.'], sep='\t')
genelist = list(set(data_gene['gene']))

cmap = pcm.get_cmap('Spectral')
line_colors = []
for i in np.arange(0, 1, 1/len(genelist)):
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
plt.savefig('%s/ahc.pdf' % save_folder)
plt.close()
nodes = R['ivl']

data_mean_feature_scale_sort = pd.DataFrame(columns=data_mean_feature_scale.columns)
for i in range(len(data_mean_feature_scale)):
    data_mean_feature_scale_sort.loc[len(data_mean_feature_scale_sort.index)] = \
        data_mean_feature_scale.iloc[int(dat.list_invert(nodes)[i])]

# heat map
plt.subplots(figsize=(8, 15))
ax1 = sns.heatmap(data_mean_feature_scale, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis')
plt.savefig('%s/heatmap.pdf' % save_folder)
plt.close()

plt.subplots(figsize=(8, 15))
ax1 = sns.heatmap(data_mean_feature_scale_sort, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis',
                  yticklabels=False)
plt.savefig('%s/heatmap_sort.pdf' % save_folder)
plt.close()

data_feature = data.copy()
all_feature = data.columns
drop_feature = [i for i in all_feature if i not in features]
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

for i in range(len(genelist)):
    if genelist[i] != 'WT':
        plt.scatter(transformed[data['gene'] == genelist[i]][0], transformed[data['gene'] == genelist[i]][1],
                    label=genelist[i], color=line_colors[i], alpha=0.5, s=1)
    else:
        plt.scatter(transformed[data['gene'] == genelist[i]][0], transformed[data['gene'] == genelist[i]][1],
                    label=genelist[i], color=line_colors[-1], alpha=0.5, s=3)
plt.legend()
plt.savefig('%s/pca.pdf' % save_folder)
plt.close()

# lda
plt.subplots(figsize=(12, 10))
lda = LDA(n_components=2)
lda_transformed = pd.DataFrame(lda.fit_transform(data_feature_scale, data['gene']))

# Plot all three series
for i in range(len(genelist)):
    if genelist[i] != 'WT':
        plt.scatter(lda_transformed[data['gene'] == genelist[i]][0], lda_transformed[data['gene'] == genelist[i]][1],
                    label=genelist[i], c=line_colors[i], alpha=0.5, s=1)
    else:
        plt.scatter(lda_transformed[data['gene'] == genelist[i]][0], lda_transformed[data['gene'] == genelist[i]][1],
                    label=genelist[i], c=line_colors[-1], alpha=0.5, s=3)

plt.legend()
plt.savefig('%s/lda.pdf' % save_folder)
plt.close()

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