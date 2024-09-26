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

batch = '5uM_48hr'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data_ctrl = data[data['treatment'] == 'DMSO'].copy().reset_index(drop=True)
print(data.head())

features = ['fov_hoechst', 'n_filtered', 'n_neg', 'n_pos', 'pos_vs_neg']

data_hp = pd.DataFrame()
for feature in features:
    data_hp['log2fc_%s' % feature] = data['log2fc_%s' % feature]
data_hp.index = data['label']

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[:50], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_0-50.pdf' % (output_dir, batch, batch))
plt.show()

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[50:100], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_50-100.pdf' % (output_dir, batch, batch))
plt.show()

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[100:150], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_100-150.pdf' % (output_dir, batch, batch))
plt.show()

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[150:200], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_150-200.pdf' % (output_dir, batch, batch))
plt.show()

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[200:250], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_200-250.pdf' % (output_dir, batch, batch))
plt.show()

fig, ax = plt.subplots(figsize=(6, 12))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp[250:], cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_250-290.pdf' % (output_dir, batch, batch))
plt.show()






