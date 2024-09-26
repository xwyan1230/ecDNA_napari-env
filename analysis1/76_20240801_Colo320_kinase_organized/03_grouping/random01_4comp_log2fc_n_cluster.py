# Fig 2b

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as shc
import utilities as uti
import matplotlib
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
mode = 'pos'  # only accepts 'neg', 'pos'
feature = 'log2fc_n_%s_mean' % mode

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
data_filter = 'qc_filter'

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    df['%s_%s' % (batch, feature)] = data[feature]
    df['%s_filter' % batch] = data[data_filter]

df_cytokinesis = pd.read_csv('%s/03_hit/02_cytokinesis.txt' % data_dir, na_values=['.'], sep='\t')
df['cytokinesis'] = df_cytokinesis['hit']

for batch in batches:
    df = df[df['%s_filter' % batch] == 1].copy().reset_index(drop=True)
    df = df.drop(['%s_filter' % batch], axis=1)
print(len(df))
df = df[df['cytokinesis'] == 0].copy().reset_index(drop=True)
print(len(df))

df_num = df.copy()
df_num = df_num.drop(['group', 'treatment', 'cytokinesis'], axis=1)

plt.figure(figsize=(20, 3))
clusters = shc.linkage(df_num, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
# plt.savefig('%s/ahc_rescale_survival_hoechst.pdf' % output_dir)
plt.show()
nodes = R['ivl']

group_seq = uti.get_seq(df['group'].tolist(), nodes)
df_group_seq = pd.DataFrame({'group': group_seq})
if not os.path.exists("%s/02_group_seq/" % output_dir):
    os.makedirs("%s/02_group_seq/" % output_dir)
df_group_seq.to_csv('%s/02_group_seq/02_survival_log2fc_n_%s.txt' % (output_dir, mode), index=False, sep='\t')

df_sort = uti.sort_df(df, group_seq)
df_ctrl = pd.DataFrame()
df_ctrl['ctrl'] = [1 if df_sort['treatment'][i] == 'DMSO' else 0 for i in range(len(df_sort))]
df_sort.index = df_sort['treatment']
df_sort = df_sort.drop(['group', 'treatment', 'cytokinesis'], axis=1)
df_sort = df_sort.astype(float)

limit = 6
fig, ax = plt.subplots(figsize=(6, 40))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(df_sort, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
# plt.savefig('%s/heatmap_%s_log2_survival_seq.pdf' % (output_dir, feature))
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

print("DONE!")