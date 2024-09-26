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
# exclude_index = [78, 85, 134, 154, 188]
exclude_index = []
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), na_values=['.'], sep='\t')

df = pd.DataFrame()
features = ['mean_n', 'log2_ratio']
for feature in features:
    for batch in batches:
        data = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]

df['survival_seq'] = survival_seq['seq']
df_drop = df.copy()
df_drop = df_drop.drop(exclude_index).reset_index(drop=True)
print(len(df_drop))
"""for batch in batches:
    df_drop = df_drop[df_drop['%s_mean_n' % batch] > 1000].reset_index(drop=True)"""
print(len(df_drop))
df_drop = df_drop.sort_values(by='survival_seq').reset_index(drop=True)
df_drop.index = df_drop['treatment']
df_drop = df_drop.drop(['index', 'treatment', 'survival_seq', 'point5uM_24hr_mean_n', 'point5uM_48hr_mean_n', '5uM_24hr_mean_n', '5uM_48hr_mean_n'], axis=1)

# heat map
fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_drop, cbar=0, linewidths=0.2, vmax=2, vmin=-2, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_ratio_survival_seq.pdf' % output_dir)
plt.show()

fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_drop[0:53], cbar=0, linewidths=0.2, vmax=2, vmin=-2, square=True, cmap='coolwarm')  #  yticklabels=False
plt.savefig('%s/heatmap_rescale_ratio_survival_seq_enlarge.pdf' % output_dir)
plt.show()