# Fig 2h

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import utilities as uti
import matplotlib
import os

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
features = ['log2fc_ratio_mean', 'log2fc_hoechst_mean', 'log2fc_n_neg_mean', 'log2fc_n_pos_mean',
            'delta_G1_mean', 'delta_G1S_mean', 'delta_S_mean', 'delta_G2M_mean', 'delta_G2MG1_mean',
            'delta_neg_G1_mean', 'delta_neg_G1S_mean', 'delta_neg_S_mean', 'delta_neg_G2M_mean', 'delta_neg_G2MG1_mean',
            'delta_pos_G1_mean', 'delta_pos_G1S_mean', 'delta_pos_S_mean', 'delta_pos_G2M_mean', 'delta_pos_G2MG1_mean',
            'log2fc_hoechst_G1_mean', 'log2fc_hoechst_G1S_mean', 'log2fc_hoechst_S_mean', 'log2fc_hoechst_G2M_mean', 'log2fc_hoechst_G2MG1_mean',
            'log2fc_hoechst_neg_G1_mean', 'log2fc_hoechst_neg_G1S_mean', 'log2fc_hoechst_neg_S_mean', 'log2fc_hoechst_neg_G2M_mean', 'log2fc_hoechst_neg_G2MG1_mean',
            'log2fc_hoechst_pos_G1_mean', 'log2fc_hoechst_pos_G1S_mean', 'log2fc_hoechst_pos_S_mean', 'log2fc_hoechst_pos_G2M_mean', 'log2fc_hoechst_pos_G2MG1_mean']
limits = [3, 0.2, 6, 6] + [0.5] * 15 + [1] * 15
filters = ['cc_filter']
filters1 = ['cytokinesis_hit']

df = pd.DataFrame()
for batch in batches:
    data = pd.read_csv('%s/01_summary/%s_grouping.txt' % (data_dir, batch), na_values=['.'], sep='\t')
    df['group'] = data['group']
    df['treatment'] = data['treatment']
    for feature in features:
        df['%s_%s' % (batch, feature)] = data[feature]
    if len(filters) != 0:
        for f in filters:
            df['%s_%s' % (batch, f)] = data[f]
    if len(filters1) != 0:
        for f in filters1:
            df['%s_%s' % (batch, f)] = data[f]
    df = df.copy()

df_flt = uti.get_df_flt(df, filters, filters1, batches)
print(len(df_flt))

df_group_seq = pd.read_csv('%s/02_group_seq/02_log2fc_ratio.txt' % data_dir, na_values=['.'], sep='\t')
group_seq = df_group_seq['group'].tolist()

df_sort = uti.sort_df(df_flt, group_seq)
df_ctrl = pd.DataFrame()
df_ctrl['ctrl'] = [1 if df_sort['treatment'][i] == 'DMSO' else 0 for i in range(len(df_sort))]

for i in range(len(features)):
    feature = features[i]
    df_feature = pd.DataFrame()
    for batch in batches:
        df_feature['%s_%s' % (batch, feature)] = df_sort['%s_%s' % (batch, feature)].tolist()
    df_feature.index = df_sort['treatment'].tolist()

    limit = limits[i]
    fig, ax = plt.subplots(figsize=(6, 40))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_feature, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/temp.pdf' % (output_dir))
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