import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%ssummary/" % master_folder

version = 2

if version == 1:
    data = pd.read_csv('%s/03_summary_log2fc_final.txt' % output_dir, na_values=['.'], sep='\t')
else:
    data = pd.read_csv('%s/03_summary_log2fc_final_v2.txt' % output_dir, na_values=['.'], sep='\t')

ctrl_AAVS = ['B6', 'C9', 'D8', 'E7', 'F2', 'G3']
ctrl_NTC = ['B10', 'C7', 'D11', 'E10', 'F9', 'G5']
ctrl_water = ['B2', 'C3', 'D4', 'E5', 'F6', 'G7']
ctrl = ctrl_AAVS + ctrl_NTC + ctrl_water
# exclude_sample = ctrl_NTC
# data = data[~data['sample'].isin(exclude_sample)].copy().reset_index(drop=True)
data_ctrl = data[data['sample'].isin(ctrl)].copy().reset_index(drop=True)

data_hp = data.copy()
data_hp.index = data['sample']
data_hp = data_hp.drop(['sample'], axis=1)

data_ctrl_hp = data_ctrl.copy()
data_ctrl_hp.index = data_ctrl['sample']
data_ctrl_hp = data_ctrl_hp.drop(['sample'], axis=1)

data_final = pd.DataFrame()

features = data_hp.columns
for feature in features:
    stdev = np.std(data_ctrl[feature])
    scaler = 0.1/stdev
    data_final[feature] = data_hp[feature] * scaler

"""data_temp = pd.DataFrame(columns = data_hp.columns)
data_temp.loc[0] = [np.nan] * len(data_hp.columns)"""

# data_final = pd.concat([data_ctrl_hp, data_temp, data_hp], axis=0)

fig, ax = plt.subplots(figsize=(7, 12))
fig.subplots_adjust(left=0.3)
ax1 = sns.heatmap(data_final, cbar=0, linewidths=2, vmax=1, vmin=-1, square=True, cmap='coolwarm', annot=False, fmt='.2f')
if version == 1:
    plt.savefig('%s/heatmap.pdf' % output_dir)
    plt.show()
else:
    plt.savefig('%s/heatmap_v2.pdf' % output_dir)
    plt.show()

data_final['sample'] = data_final.index
if version == 1:
    data_final.to_csv('%s/04_summary_normalized.txt' % output_dir, index=False, sep='\t')
else:
    data_final.to_csv('%s/04_summary_normalized_v2.txt' % output_dir, index=False, sep='\t')