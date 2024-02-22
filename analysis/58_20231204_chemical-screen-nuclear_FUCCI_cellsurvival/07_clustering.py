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

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231204_analysis_FUCCI-screen/"
df = pd.read_excel('%slog2FC_summary.xlsx' % master_folder)

ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']
sample = ['B8', 'C9', 'E6', 'G2', 'B7', 'B10', 'C7',
          'D5', 'C4', 'B5',
          'E4', 'B9', 'G7', 'G6', 'E10', 'D4', 'G10', 'B4', 'C11',
          'D9', 'E2', 'E9', 'F11', 'G5', 'G11']
# rank by DM_24hr
sample = ['G2', 'B7', 'C9', 'B10', 'E6', 'C7', 'B8', 'D5', 'C4', 'B5', 'D4', 'G7', 'E10', 'E4', 'G10', 'B9', 'G6', 'B4', 'C11', 'F5', 'D9', 'E5', 'B6', 'B3', 'E7', 'B2', 'D8', 'F9', 'B11', 'G11', 'C8', 'D10', 'G5', 'E11', 'C6', 'G9', 'G8', 'D2', 'F8', 'E9', 'D11', 'E3', 'E8', 'F4', 'C2', 'C5', 'D7', 'F7', 'G4', 'G3', 'E2', 'F11', 'F6', 'F2', 'D3']
# rank by DM_24hr + HSR_24hr
# sample = ['C9', 'G2', 'B8', 'B7', 'E6', 'B10', 'C4', 'C7', 'B5', 'E4', 'B9', 'D5', 'G7', 'D4', 'E10', 'G6', 'G10', 'B4', 'D8', 'C11', 'F5', 'D9', 'E7', 'B6', 'B2', 'E5', 'B3', 'B11', 'F9', 'F4', 'C6', 'G11', 'D7', 'C8', 'D11', 'E9', 'D10', 'F7', 'F10', 'D2', 'F8', 'G3', 'G9', 'E8',  'C5', 'F6', 'G5', 'G4', 'E11', 'C2', 'G8', 'E2', 'D3', 'F2', 'E3']
# DM only
sample = ['B10', 'C7', 'D5', 'D4', 'G10', 'B4', 'C11', 'F5', 'E5', 'B3', 'F9', 'B11']
# both
sample = ['B8', 'E4', 'D8']
# DM sensitive
sample = ['G2', 'B7', 'C9', 'E6', 'C4', 'B5', 'G7', 'E10', 'B9', 'G6', 'D9', 'B6', 'E7', 'B2']

samples = ctrl + sample

df_sample = pd.DataFrame()
for i in range(len(samples)):
    df_temp = df[df['sample'] == samples[i]].copy().reset_index(drop=True)
    df_sample = pd.concat([df_sample, df_temp], axis=0)
df_sample.index = df_sample['sample']
df_sample = df_sample.drop('sample', axis=1)

fig, ax = plt.subplots(figsize=(9, 14))
fig.subplots_adjust(left=0.3)
ax1 = sns.heatmap(df_sample, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f') # annot=True
plt.savefig('%s/heatmap_DMsens_rankby_DM_24hr.pdf' % master_folder)
plt.show()