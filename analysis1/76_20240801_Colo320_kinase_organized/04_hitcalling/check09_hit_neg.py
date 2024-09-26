# Fig 2e

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as shc
import utilities as uti
import matplotlib
import numpy as np
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
mode = 'pos'

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
feature = 'log2fc_n_%s_mean' % mode
target_df = pd.read_csv('%s/targets.txt' % data_dir1, na_values=['.'], sep='\t')
n_target_limit = 5

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_hitcalling.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    data1 = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    df['%s_%s_hit' % (batch, feature)] = data[feature]
    df['%s_%s' % (batch, feature)] = data1[feature]

for batch in batches:
    df = df[df['%s_%s_hit' % (batch, feature)] != -100].copy().reset_index(drop=True)
print(len(df))

df['hit_sum'] = np.sum(df[['%s_%s_hit' % (batch, feature) for batch in batches]], axis=1)
df['sum'] = np.sum(df[['%s_%s' % (batch, feature) for batch in batches]], axis=1)
sum_lst = []
for i in range(len(df)):
    temp_lst = [df['point5uM_24hr_%s' % feature][i], df['point5uM_48hr_%s' % feature][i], df['5uM_24hr_%s' % feature][i], df['5uM_48hr_%s' % feature][i]]
    temp_lst.sort(reverse=True)
    print(temp_lst)
    sum_lst.append(temp_lst[0] + temp_lst[1] + temp_lst[2] + temp_lst[3])
df['sum'] = sum_lst

df_num = df.copy()
df_num = df_num[df_num['hit_sum']<=-1]
df_num = df_num.sort_values(by=['hit_sum', 'sum'], ascending=[True, True])

print(len(df_num[df_num['hit_sum'] == -4]))
print(len(df_num[df_num['hit_sum'] == -3]))
print(len(df_num[df_num['hit_sum'] == -2]))
print(len(df_num[df_num['hit_sum'] == -1]))
print(df_num)

df_num1 = pd.DataFrame()
for batch in batches:
    df_num1['%s_%s' % (batch, feature)] = df_num['%s_%s' % (batch, feature)]
df_num1.index = df_num['treatment']
limit = 3
fig, ax = plt.subplots(figsize=(6, 20))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_num1, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
# plt.savefig('%s/heatmap_ratio.pdf' % (output_dir))
plt.show()

"""hits = df_num[df_num['hit_sum'] >= 2]['treatment'].tolist()
print(len(hits))

df_sort = df.copy()
df_ctrl = pd.DataFrame()
df_sort = df_sort.sort_values(by=['sum', 'hit_sum'], ascending=[True, True]).reset_index(drop=True)
df_ctrl['ctrl'] = [1 if df_sort['treatment'][i] == 'DMSO' else 0 for i in range(len(df_sort))]
df_num = pd.DataFrame()
for batch in batches:
    df_num['%s_%s' % (batch, feature)] = df_sort['%s_%s' % (batch, feature)]
df_num.index = df_sort['treatment']

df_group_seq = pd.DataFrame({'group': df_sort['group'].tolist()})
df_group_seq.to_csv('%s/02_group_seq/02_log2fc_ratio.txt' % (output_dir), index=False, sep='\t')

limit = 3
fig, ax = plt.subplots(figsize=(6, 80))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_num, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
# plt.savefig('%s/heatmap_ratio.pdf' % (output_dir))
plt.show()

norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), "darkblue"],
          [norm(-0.6), "#d0d2d3"],
          [norm( 0.6), "#d0d2d3"],
          [norm( 1.0), rainboo_colors[6]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_ctrl, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap, yticklabels=False)  #  yticklabels=False
# plt.savefig('%s/heatmap_rescale_survival_hoechst_DMSO.pdf' % output_dir)
plt.show()


data, data_n, data_percentage = uti.get_gene_stats(hits, target_df, n_target_limit)
print(data)

norm = matplotlib.colors.Normalize(-1, 1)
colors = [[norm(-1.0), "#d0d2d3"],
          [norm(0), "#d0d2d3"],
          [norm(1.0), rainboo_colors[0]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

limit = 4
fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_n, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True,
                  cmap=cmap)  # yticklabels=False
# plt.savefig('%s/heatmap_%s_log2_survival_seq_hit.pdf' % (output_dir, feature))
plt.show()

fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_percentage, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True,
                  cmap=cmap)  # yticklabels=False
# plt.savefig('%s/heatmap_%s_log2_survival_seq_hit.pdf' % (output_dir, feature))
plt.show()
"""