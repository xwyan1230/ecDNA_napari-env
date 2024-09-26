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

search = 'DMSO'
data_flt = data[data['target'].str.contains(search)].copy().reset_index(drop=True)

features = ['fov_hoechst', 'n_filtered', 'n_neg', 'n_pos','pos_vs_neg']

data_hp = pd.DataFrame()
for feature in features:
    data_hp['log2fc_%s' % feature] = data_flt['log2fc_%s' % feature]
data_hp.index = data_flt['label']

fig, ax = plt.subplots(figsize=(5, len(data_flt) * 0.3 + 1))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_%s.pdf' % (output_dir, batch, batch, search))
plt.show()







