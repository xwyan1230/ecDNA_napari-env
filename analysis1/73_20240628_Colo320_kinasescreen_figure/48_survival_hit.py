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

data = pd.read_csv('%s/common_hit_percentage.txt' % (output_dir), na_values=['.'], sep='\t')
# data = data[data['percentage']>0.5]
data = data.sort_values(by=['n', 'percentage'], ascending=[False, False]).reset_index(drop=True)
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
print(data)

norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1), "#d0d2d3"], # "darkblue"
          [norm(0), "#d0d2d3"],
          [norm( 0), "#d0d2d3"],
          [norm( 1), rainboo_colors[0]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

df = data.copy()
df.index = df['gene']
df = df.drop(['gene', 'percentage'], axis=1)
fig, ax = plt.subplots(figsize=(6, 8))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df, linewidths=0.2, vmax=4, vmin=-4, square=True, cmap=cmap)  #  yticklabels=False
plt.savefig('%s/heatmap_common_survival_hit_gene_n.pdf' % output_dir)
plt.show()

df = data.copy()
df.index = df['gene']
df = df.drop(['gene', 'n'], axis=1)
fig, ax = plt.subplots(figsize=(6, 8))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap)  #  yticklabels=False
plt.savefig('%s/heatmap_common_survival_hit_gene_percentage.pdf' % output_dir)
plt.show()