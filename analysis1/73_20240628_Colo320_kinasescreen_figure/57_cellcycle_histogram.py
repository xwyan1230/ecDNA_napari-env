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
    # data['delta_sum'] = np.abs(data['mean_delta_n_G1']) + np.abs(data['mean_delta_n_G1S']) + np.abs(data['mean_delta_n_S']) + np.abs(data['mean_delta_n_G2M']) + np.abs(data['mean_delta_n_G2MG1'])
    data1 = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
    df_temp = pd.DataFrame()
    df_temp['index'] = data.index
    df_temp['batch'] = [batch] * len(df_temp)
    # df_temp['delta_sum'] = data['delta_sum']
    df_temp['G1'] = data['mean_per_n_G1']
    df_temp['G1S'] = data['mean_per_n_G1S']
    df_temp['S'] = data['mean_per_n_S']
    df_temp['G2M'] = data['mean_per_n_G2M']
    df_temp['G2MG1'] = data['mean_per_n_G2MG1']
    df_temp['treatment'] = data['treatment']
    df_temp['log2_mean_hoechst'] = data1['log2_mean_hoechst']
    df = pd.concat([df, df_temp], axis=0).reset_index(drop=True)

interval = [-0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0]
group = []
pointer = 0
for i in range(len(df)):
    for j in range(len(interval)):
        if (df['log2_mean_hoechst'][i] < interval[j]) & (pointer == 0):
            group.append(j)
            pointer = 1
    if pointer == 0:
        group.append(8)
    else:
        pointer = 0
df['group'] = group
df_drop = df.copy()
df_drop = df_drop[~df_drop['index'].isin(exclude_index)].reset_index(drop=True)
print(len(df_drop))

cc = ['G1', 'G1S', 'S', 'G2M', 'G2MG1']

"""df_cc = pd.DataFrame()
df_cc['group'] = np.arange(0, 9, 1).tolist() * 5
df_cc['cellcycle'] = ['G1'] * 9 + ['G1S'] * 9 + ['S'] * 9 + ['G2M'] * 9 + ['G2MG1'] * 9
avg_lst = []
for j in range(5):
    for i in range(9):
        avg_lst.append(np.mean(df[df['group'] == i][cc[j]]))
df_cc['percentage'] = avg_lst"""

df_cc = pd.DataFrame()
df_cc['group'] = np.arange(0, 9, 1)
for j in range(5):
    avg_lst = []
    for i in range(9):
        avg_lst.append(np.mean(df[df['group'] == i][cc[j]]))
    df_cc[cc[j]] = avg_lst

print(df_cc)

# fig, ax = plt.subplots(figsize=(9, 7))
"""# fig.subplots_adjust(right=0.8)
sns.barplot(data=df_cc, x='group', y='percentage', hue='cellcycle', stacked=True)  ##4d4e4e
# plt.savefig('%s/cell-cycle_vs_survival_update1.pdf' % (output_dir))
plt.show()"""

colors = ['#bc4d4a', '#3d5a80', '#669daa', '#dd933c', '#d2b48c'] ## Colors for the bars
df_cc.set_index('group').plot(kind='bar', stacked=True, color=colors) ## Plot
plt.ticklabel_format(style='plain', useOffset=False, axis='y') ## No offset
plt.gca().set_ylabel("Total $'s")
plt.savefig('%s/cell-cycle_group_by_survival.pdf' % (output_dir))
plt.show()