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
hit = 'survival_all'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]

hits = pd.read_csv('%s/hit.txt' % master_folder, na_values=['.'], sep='\t')

if hit == 'survival_all':
    features = ['log2fc_fov_hoechst', 'log2fc_n_filtered']
    data_flt = data[(data[features[0]] <= -1) | (data[features[0]] >= 1)].copy().reset_index(drop=True)
    """xs = data_flt.columns
    print(xs)
    import collections
    print([item for item, count in collections.Counter(xs).items() if count > 1])"""
    for i in range(len(features)-1):
        data_temp = data[(data[features[i+1]] <= -1) | (data[features[i+1]] >= 1)].copy().reset_index(drop=True)
        data_flt = data_flt.merge(data_temp, how='inner').reset_index(drop=True)
else:
    if hit == 'survival':
        features = ['log2fc_fov_hoechst']
    elif hit == 'survival_n':
        features = ['log2fc_n_filtered']
    elif hit == 'cytokinesis':
        features = ['log2fc_cytokinesis']
    else:
        features = ['per_G1', 'per_G1S', 'per_S', 'per_G2M']

    data_flt = pd.DataFrame()
    for feature in features:
        data_temp = data[(data[feature] <= -1) | (data[feature] >= 1)]
        data_flt = pd.concat([data_flt, data_temp], axis=0).drop_duplicates().reset_index(drop=True)

if hit in ['survival', 'survival_all']:
    data_flt = data_flt.sort_values(by=['log2fc_fov_hoechst']).reset_index(drop=True)
elif hit == 'survival_n':
    data_flt = data_flt.sort_values(by=['log2fc_n_filtered']).reset_index(drop=True)
elif hit == 'cytokinesis':
    data_flt1 = data_flt[data_flt['mean_n_filtered'] >= 200].copy().reset_index(drop=True)
    data_flt = data_flt1.sort_values(by=['log2fc_cytokinesis'], ascending=False).reset_index(drop=True)
else:
    data_flt = data_flt.sort_values(by=['cc_score']).reset_index(drop=True)

data_flt_na = data_flt.copy()
features = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S', 'log2fc_per_G2M', 'log2fc_cytokinesis']
for feature in features:
    data_flt_na.loc[(data_flt_na['rep'] == 2) & (data_flt_na['mean_n_filtered'] < 200), feature] = np.nan

features = ['log2fc_fov_hoechst', 'log2fc_n_filtered', 'cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S',
            'log2fc_per_G2M', 'log2fc_cytokinesis']

data_hp = pd.DataFrame()
data_hp_na = pd.DataFrame()
for feature in features:
    data_hp[feature] = data_flt[feature]
    data_hp_na[feature] = data_flt_na[feature]
data_hp.index = data_flt['label']
data_hp_na.index = data_flt_na['label']

hits.loc[len(hits.index)] = [batch, hit, len(data_flt['label'].tolist()), len(data_flt['label'].tolist())/240, str(['%s*' % data_flt['label'].tolist()[i] for i in range(len(data_flt))])]
print(hits.head())
hits = hits.drop_duplicates().reset_index(drop=True)
hits.to_csv('%s/hit.txt' % master_folder, index=False, sep='\t')

fig, ax = plt.subplots(figsize=(7, len(data_flt) * 0.3 + 1))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_global_%s.pdf' % (output_dir, batch, batch, hit))
plt.show()

fig, ax = plt.subplots(figsize=(7, len(data_flt) * 0.3 + 1))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_hp_na, cbar=0, linewidths=2, vmax=4, vmin=-4, square=True, cmap='coolwarm', annot=False, fmt='.2f')
plt.savefig('%s/%s/%s_heatmap_global_%s_na.pdf' % (output_dir, batch, batch, hit))
plt.show()


