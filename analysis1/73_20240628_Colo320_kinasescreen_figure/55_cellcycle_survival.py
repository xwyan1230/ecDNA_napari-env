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
from scipy.stats import pearsonr
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
exclude_index = [154, 188]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77',
                  '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/%s/%s_cellcycle_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    data['delta_sum'] = np.abs(data['mean_delta_n_G1']) + np.abs(data['mean_delta_n_G1S']) + np.abs(data['mean_delta_n_S']) + np.abs(data['mean_delta_n_G2M']) + np.abs(data['mean_delta_n_G2MG1'])
    data1 = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    df_temp = pd.DataFrame()
    df_temp['index'] = data.index
    df_temp['batch'] = [batch] * len(df_temp)
    df_temp['delta_sum'] = data['delta_sum']
    df_temp['treatment'] = data['treatment']
    df_temp['log2_mean_hoechst'] = data1['log2_mean_hoechst']
    df = pd.concat([df, df_temp], axis=0).reset_index(drop=True)

df_drop = df.copy()
df_drop = df_drop[~df_drop['index'].isin(exclude_index)].reset_index(drop=True)
print(len(df_drop))
fig, ax = plt.subplots(figsize=(9, 7))
# fig.subplots_adjust(right=0.8)
sns.scatterplot(data=df_drop, x='log2_mean_hoechst', y='delta_sum', s=8, alpha=1, color='#cccccc')  ##4d4e4e
sns.scatterplot(data=df_drop[df_drop['treatment'] == 'DMSO'], x='log2_mean_hoechst', y='delta_sum', s=8, alpha=1, color=rainboo_colors[6])
plt.savefig('%s/cell-cycle_vs_survival_update1.pdf' % (output_dir))
plt.show()

corr, _ = pearsonr(df_drop['log2_mean_hoechst'], df_drop['delta_sum'])
print('Pearsons correlation: %.3f' % corr)