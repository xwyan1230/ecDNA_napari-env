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
print(df.head())
ctrl = ['C3', 'C10', 'D6', 'F3', 'F10']

"""fig, ax = plt.subplots(figsize=(10, 10))
# fig.subplots_adjust(right=0.8)
# line_colors = [(255 / 255, 222 / 255, 173 / 255)] * 6
# sns.set_palette(sns.color_palette(line_colors))
sns.scatterplot(data=df, x='HSR_24hr_per_survival_log2FC', y='DM_24hr_per_survival_log2FC')
plt.legend(loc=(1.04, 0))
plt.xlim([-9, 1])
plt.ylim([-9, 1])
# plt.savefig('%s/%s_growth.pdf' % (master_folder, wells[w]))
plt.show()"""

"""plt.subplots(figsize=(8, 8))
# plt.scatter(df['HSR_2hr_per_survival_log2FC'], df['DM_2hr_per_survival_log2FC'], color=(255/255, 215/255, 0/255), s=20, alpha=0.5)
# plt.scatter(df['HSR_6hr_per_survival_log2FC'], df['DM_6hr_per_survival_log2FC'], color=(0/255, 191/255, 255/255), s=20, alpha=0.5)
plt.scatter(df['HSR_24hr_per_survival_log2FC'], df['DM_24hr_per_survival_log2FC'], color=(154/255, 205/255, 50/255), s=20, alpha=0.5)
df_ctrl = df[df['sample'].isin(ctrl)].copy().reset_index(drop=True)
# plt.scatter(df_ctrl['HSR_2hr_per_survival_log2FC'], df_ctrl['DM_2hr_per_survival_log2FC'], color=(220/255, 20/255, 60/255), s=20, alpha=0.5)
# plt.scatter(df_ctrl['HSR_6hr_per_survival_log2FC'], df_ctrl['DM_6hr_per_survival_log2FC'], color=(220/255, 20/255, 60/255), s=20, alpha=0.5)
plt.scatter(df_ctrl['HSR_24hr_per_survival_log2FC'], df_ctrl['DM_24hr_per_survival_log2FC'], color=(220/255, 20/255, 60/255), s=20, alpha=0.5)
# df_temp_6hr = df[df['DM_6hr_per_survival_log2FC'] < -0.3].copy().reset_index(drop=True)
# for i in range(len(df_temp_6hr)):
#     plt.text(x=df_temp_6hr['HSR_6hr_per_survival_log2FC'][i]+0.005, y=df_temp_6hr['DM_6hr_per_survival_log2FC'][i]+0.005, s=df_temp_6hr['sample'][i], size=5, color=(0/255, 191/255, 255/255))
df_temp_24hr = df[df['DM_24hr_per_survival_log2FC'] < -0.3].copy().reset_index(drop=True)
for i in range(len(df_temp_24hr)):
    plt.text(x=df_temp_24hr['HSR_24hr_per_survival_log2FC'][i]+0.005, y=df_temp_24hr['DM_24hr_per_survival_log2FC'][i]+0.005, s=df_temp_24hr['sample'][i], size=5, color=(154/255, 205/255, 50/255))
plt.xlim([-9, 1])
plt.ylim([-9, 1])
plt.xlabel('HSR')
plt.ylabel('DM')
plt.savefig('%s/growth_24hr.pdf' % master_folder)
plt.show()"""

df['new'] = df['DM_24hr_per_survival_log2FC'] + df['HSR_24hr_per_survival_log2FC']
df_temp = df.sort_values(by='new')
print(df_temp['sample'].tolist())